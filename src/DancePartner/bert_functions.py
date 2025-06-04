import os 
import pandas as pd
from transformers import BertTokenizer, TrainingArguments, Trainer
from torch import nn
import numpy as np
import re

from sklearn.preprocessing import LabelEncoder
from transformers import BertTokenizer, TrainingArguments, Trainer

import torch
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, Dropout, Linear, MSELoss, Tanh
from transformers import BertModel, BertPreTrainedModel
from transformers.modeling_outputs import SequenceClassifierOutput

from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, precision_recall_fscore_support

def __make_bert_ready(df: pd.DataFrame, segment_col_name: str):
    """
    Function to prepare a dataframe to be inputted into the BERT model

    Parameters
    ----------
    df 
        The dataframe of setnences. Should be a result of `find_names_in_papers`
    
    segment_col_name
        The name of the column representing the chunk of text containing the pair of biomolecules.
    
    Returns
    -------
        A Pandas.DataFrame that has been cleaned
    """

    # Add required columns to data frame. BERT cannot handle more than 250 characters 
    df['Term1'] = "@TERM$1"
    df['Term2'] = "@TERM$2"
    df_clean = df[(df[segment_col_name].str.len() < 250) & (df[segment_col_name].str.len() > 0)] #Check that text is btwn 1 and 250 characters
    
    # Replace first occurences of protein terms with MASK variables
    df_clean['Sentence'] = df_clean.apply(lambda row: re.sub(str(row["term_1"]), " @TERM$1 ", row[segment_col_name], count=1), axis=1)
    df_clean['Sentence'] = df_clean.apply(lambda row: re.sub(str(row["term_2"]), " @TERM$2 ", row.Sentence, count=1), axis=1)
    df_clean['Guess'] = 0
    
    # Check that '@TERM$1' and '@TERM$2' only appear once in each sentence segment
    df_checked = df_clean[df_clean.Sentence.str.count("@TERM\$1") == 1]
    df_checked = df_checked[df_checked.Sentence.str.count("@TERM\$2") == 1]
    df_checked = df_checked.drop(columns=[segment_col_name])
    return(df_checked)

def __preprocess_data(df: pd.DataFrame, tokenizer: BertTokenizer, x_col: str, y_col: str, e1_col: str, e2_col: str):
    """
    Function to preprocess data for BERT

    Parameters
    ----------
    df
        The dataframe from __make_bert_ready
    
    tokenizer
        A tokenizer object from the transformers package

    x_col
        The sentence column name

    y_col
        The label column name

    e1_col
        The first entity column name

    e2_col
        The second entity column name

    Returns
    -------
        A SrcDataset
    """
    # 2-Masted-senteces input format
    x1 = df.apply(lambda x: x[x_col].replace(x[e1_col], "[MASK]"), axis=1).tolist()
    x2 = df.apply(lambda x: x[x_col].replace(x[e2_col], "[MASK]"), axis=1).tolist()

    label_encoder = LabelEncoder()
    y = torch.tensor(label_encoder.fit_transform(df[y_col]), dtype=torch.long)

    tokenized_x = tokenizer(x1, x2, return_tensors="pt", padding='max_length', truncation=True, max_length=512)

    dataset = SrcDataset(tokenized_x, y)
    return dataset

class SrcDataset(torch.utils.data.Dataset):
    """
    Dataset class for BertSRC
    """
    def __init__(self, encodings, labels):
        self.encodings = encodings
        self.labels = labels

    def __getitem__(self, idx):
        item = {key: val[idx] for key, val in self.encodings.items()}
        item["labels"] = self.labels[idx]
        return item

    def __len__(self):
        return len(self.labels)

class BertSrcClassifier(BertPreTrainedModel):
    """
    Dataset class for BertSRC Classifier
    """
    def __init__(self, config, mask_token_id: int, num_token_layer: int = 2):
        super().__init__(config)
        self.mask_token_id = mask_token_id
        self.n_output_layer = num_token_layer
        self.num_labels = config.num_labels
        self.config = config

        self.bert = BertModel(config)
        self.dense = Linear(config.hidden_size * num_token_layer, config.hidden_size)
        self.activation = Tanh()
        self.dropout = Dropout(p = 0.5) # Originally 0.1
        self.classifier = Linear(config.hidden_size, config.num_labels)

        self.post_init()

    def forward(
        self,
        input_ids: torch.Tensor = None,
        attention_mask: torch.Tensor = None,
        token_type_ids: torch.Tensor = None,
        position_ids: torch.Tensor = None,
        head_mask: torch.Tensor = None,
        inputs_embeds: torch.Tensor = None,
        labels: torch.Tensor = None,
        output_attentions: bool = None,
        output_hidden_states: bool = None,
        return_dict: bool = None,
    ):
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.bert(
            input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        assert 1 <= self.n_output_layer <= 3
        if self.n_output_layer == 1:
            output = outputs[0][0]
        else:
            check = input_ids == self.mask_token_id
            if self.n_output_layer == 3:
                check[:, 0] = True
            output = torch.reshape(
                outputs[0][check], (-1, self.n_output_layer * self.config.hidden_size)
            )

        output = self.dense(output)
        output = self.activation(output)
        output = self.dropout(output)
        logits = self.classifier(output)

        loss = None
        if labels is not None:
            if self.config.problem_type is None:
                if self.num_labels == 1:
                    self.config.problem_type = "regression"
                elif self.num_labels > 1 and (
                    labels.dtype == torch.long or labels.dtype == torch.int
                ):
                    self.config.problem_type = "single_label_classification"
                else:
                    self.config.problem_type = "multi_label_classification"

            if self.config.problem_type == "regression":
                loss_fct = MSELoss()
                if self.num_labels == 1:
                    loss = loss_fct(logits.squeeze(), labels.squeeze())
                else:
                    loss = loss_fct(logits, labels)
            elif self.config.problem_type == "single_label_classification":
                loss_fct = CrossEntropyLoss()
                loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))
            elif self.config.problem_type == "multi_label_classification":
                loss_fct = BCEWithLogitsLoss()
                loss = loss_fct(logits, labels)
        if not return_dict:
            output = (logits,) + outputs[2:]
            return ((loss,) + output) if loss is not None else output

        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )

def run_bert(input_path: str, model_path: str, output_directory: str, segment_col_name: str, **kwargz):
    """
    Function to prepare a dataframe to be inputted into the BERT model

    Parameters
    ----------
    input_path
        A path to the CSV file to run the model on. Should be a result of `ppi.find_terms_in_papers`
    
    model_path
        A path to the folder containing the BERT model. Put the model in a folder within this directory called biobert. 
        Find the model here: https://huggingface.co/david-degnan/BioBERT-RE/tree/main
    
    output_directory
        A path where to write the results to
    
    segment_col_name
        The name of the column representing the chunk of text containing the pair of biomolecules.
    
    **kwargz
        Any additional arguments to pass to `TrainingArguments`.
    
    Returns
    -------
        Writes a csv file containing the results of the model.
    """

    CLASSES = ["0", "1"]

    tokenizer = BertTokenizer.from_pretrained(
        "dmis-lab/biobert-base-cased-v1.2"
    )

    model = BertSrcClassifier.from_pretrained(
        model_path,
        num_labels=len(CLASSES),
        mask_token_id=tokenizer.mask_token_id,
    )

    test = pd.read_csv(input_path)
    test = __make_bert_ready(test, segment_col_name)
    test_dataset = __preprocess_data(test, tokenizer, x_col = "Sentence", y_col = "Guess", e1_col = "Term1", e2_col = "Term2")

    training_args = TrainingArguments(
        output_dir="./checkpoints",
        logging_dir="./logs",
        **kwargz
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        eval_dataset=test_dataset,
    )

    trainer_out = trainer.predict(test_dataset)

    probs = nn.functional.softmax(torch.from_numpy(trainer_out.predictions), dim = -1)
    test[["True Negative", "True Positive"]] = pd.DataFrame(probs)
    test = test.drop("Guess", axis = 1)
    test.to_csv(os.path.join(output_directory, "bert_results.txt"), sep = '\t', index = False, header = True)

    return(None)