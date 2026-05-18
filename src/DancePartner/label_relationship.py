import pandas as pd

# Define each term and its relationship type (positive, negative, unknown) to be used in labeling relationships in the BERT output.
term_relationships = {
    "abate": "negative", "abated": "negative", "abating": "negative", "abatement": "negative", "abolish": "negative",
    "abolished": "negative", "abolishing": "negative", "abolition": "negative", "acetylate": "unknown", "acetylated": "unknown",
    "acetylating": "unknown", "acetylation": "unknown", "acrylate": "unknown", "acrylated": "unknown", "acrylating": "unknown",
    "acrylation": "unknown", "activate": "positive", "activated": "positive", "activating": "positive", "activation": "positive",
    "acylate": "unknown", "acylated": "unknown", "acylating": "unknown", "acylation": "unknown", "adhere": "unknown",
    "adhered": "unknown", "adhering": "unknown", "adhesion": "unknown", "affix": "unknown", "affixed": "unknown",
    "affixing": "unknown", "affixion": "unknown", "aggregate": "unknown", "aggregated": "unknown", "aggregating": "unknown",
    "aggregation": "unknown", "align": "unknown", "aligned": "unknown", "aligning": "unknown", "alignment": "unknown",
    "alkylate": "unknown", "alkylated": "unknown", "alkylating": "unknown", "alkylation": "unknown", "anabolize": "positive",
    "anabolized": "positive", "anabolizing": "positive", "anabolism": "positive", "annex": "unknown", "annexed": "unknown",
    "annexing": "unknown", "annexation": "unknown", "append": "unknown", "appended": "unknown", "appending": "unknown",
    "appendage": "unknown", "assemble": "unknown", "assembled": "unknown", "assembling": "unknown", "assembly": "unknown",
    "associate": "unknown", "associated": "unknown", "associating": "unknown", "association": "unknown", "attach": "unknown",
    "attached": "unknown", "attaching": "unknown", "attachment": "unknown", "attenuate": "negative", "attenuated": "negative",
    "attenuating": "negative", "attenuation": "negative", "bind": "unknown", "binding": "unknown", "block": "negative",
    "blocked": "negative", "blocking": "negative", "blockage": "negative", "bound": "unknown", "bridge": "unknown",
    "bridged": "unknown", "bridging": "unknown", "butylate": "unknown", "butylated": "unknown", "butylating": "unknown",
    "butylation": "unknown", "carboxylate": "unknown", "carboxylated": "unknown", "carboxylating": "unknown", "carboxylation": "unknown",
    "catabolize": "negative", "catabolized": "negative", "catabolizing": "negative", "catabolism": "negative", "catalyze": "unknown",
    "catalysis": "unknown", "catalyzed": "unknown", "catalyzing": "unknown", "change": "unknown", "changed": "unknown",
    "changeover": "unknown", "changing": "unknown", "chain": "unknown", "chained": "unknown", "chaining": "unknown",
    "cleave": "negative", "cleavage": "negative", "cleaved": "negative", "cleaving": "negative", "cluster": "unknown",
    "clustered": "unknown", "clustering": "unknown", "cohere": "unknown", "cohered": "unknown", "cohering": "unknown",
    "cohesion": "unknown", "combine": "unknown", "combined": "unknown", "combining": "unknown", "combination": "unknown",
    "complex": "unknown", "complexation": "unknown", "complexed": "unknown", "complexing": "unknown", "confine": "negative",
    "confined": "negative", "confining": "negative", "confinement": "negative", "connect": "unknown", "connected": "unknown",
    "connecting": "unknown", "connection": "unknown", "constrain": "negative", "constrained": "negative", "constraining": "negative",
    "constraint": "negative", "constrict": "negative", "constricted": "negative", "constricting": "negative", "constriction": "negative",
    "couple": "unknown", "coupled": "unknown", "coupling": "unknown", "create": "positive", "created": "positive",
    "creating": "positive", "creation": "positive", "crylate": "unknown", "crylated": "unknown", "crylating": "unknown",
    "crylation": "unknown", "decrease": "negative", "decreased": "negative", "decreasing": "negative", "decrement": "negative",
    "detach": "negative", "detached": "negative", "detaching": "negative", "detachment": "negative", "detain": "unknown",
    "detained": "unknown", "detaining": "unknown", "detention": "unknown", "deter": "negative", "deterred": "negative",
    "deterrence": "negative", "deterring": "negative", "dilute": "negative", "diluted": "negative", "diluting": "negative",
    "dilution": "negative", "dimerization": "unknown", "dimerize": "unknown", "dimerized": "unknown", "dimerizing": "unknown",
    "diminish": "negative", "diminished": "negative", "diminishment": "negative", "diminishing": "negative", "disassemble": "negative",
    "disassembled": "negative", "disassembling": "negative", "disassembly": "negative", "dissipate": "negative", "dissipated": "negative",
    "dissipating": "negative", "dissipation": "negative", "elevate": "positive", "elevated": "positive", "elevating": "positive",
    "elevation": "positive", "eliminate": "negative", "eliminated": "negative", "eliminating": "negative", "elimination": "negative",
    "enhance": "positive", "enhanced": "positive", "enhancing": "positive", "enhancement": "positive", "ethylate": "unknown",
    "ethylated": "unknown", "ethylating": "unknown", "ethylation": "unknown", "extenuate": "negative", "extenuated": "negative",
    "extenuating": "negative", "extenuation": "negative", "facilitate": "positive", "facilitated": "positive", "facilitating": "positive",
    "facilitation": "positive", "fasten": "unknown", "fastened": "unknown", "fastening": "unknown", "free": "unknown",
    "freedom": "unknown", "freed": "unknown", "freeing": "unknown", "fuse": "unknown", "fused": "unknown",
    "fusing": "unknown", "fusion": "unknown", "generate": "positive", "generated": "positive", "generating": "positive",
    "generation": "positive", "glycosylate": "unknown", "glycosylated": "unknown", "glycosylating": "unknown", "glycosylation": "unknown",
    "group": "unknown", "grouped": "unknown", "grouping": "unknown", "hemolysis": "negative", "hemolyze": "negative",
    "hemolyzed": "negative", "hemolyzing": "negative", "hinder": "negative", "hindered": "negative", "hindering": "negative",
    "hindrance": "negative", "hydrolysis": "negative", "hydrolyze": "negative", "hydrolyzed": "negative", "hydrolyzing": "negative",
    "impair": "negative", "impaired": "negative", "impairing": "negative", "impairment": "negative", "impedance": "negative",
    "impede": "negative", "impeded": "negative", "impeding": "negative", "increase": "positive", "increased": "positive",
    "increasing": "positive", "increment": "positive", "induce": "positive", "induced": "positive", "inducing": "positive",
    "induction": "positive", "inhibit": "negative", "inhibited": "negative", "inhibiting": "negative", "inhibition": "negative",
    "interact": "unknown", "interacted": "unknown", "interacting": "unknown", "interaction": "unknown", "intercept": "negative",
    "intercepted": "negative", "intercepting": "negative", "interception": "negative", "interfere": "negative", "interference": "negative",
    "interfered": "negative", "interfering": "negative", "intricate": "unknown", "intricated": "unknown", "intricating": "unknown",
    "intrication": "unknown", "join": "unknown", "joined": "unknown", "joining": "unknown", "liberate": "unknown",
    "liberated": "unknown", "liberating": "unknown", "liberation": "unknown", "ligate": "unknown", "ligated": "unknown",
    "ligating": "unknown", "ligation": "unknown", "link": "unknown", "linkage": "unknown", "linked": "unknown",
    "linking": "unknown", "loosen": "negative", "loosened": "negative", "loosening": "negative", "metabolize": "unknown",
    "metabolized": "unknown", "metabolizing": "unknown", "metabolism": "unknown", "methylate": "unknown", "methylated": "unknown",
    "methylating": "unknown", "methylation": "unknown", "moderate": "unknown", "moderated": "unknown", "moderating": "unknown",
    "moderation": "unknown", "modulate": "unknown", "modulated": "unknown", "modulating": "unknown", "modulation": "unknown",
    "neutralization": "unknown", "neutralize": "unknown", "neutralized": "unknown", "neutralizing": "unknown", "obstruct": "negative",
    "obstructed": "negative", "obstructing": "negative", "obstruction": "negative", "occlude": "negative", "occluded": "negative",
    "occluding": "negative", "occlusion": "negative", "oligomerization": "unknown", "oligomerize": "unknown", "oligomerized": "unknown",
    "oligomerizing": "unknown", "organization": "unknown", "organize": "unknown", "organized": "unknown", "organizing": "unknown",
    "osmolysis": "negative", "osmolyze": "negative", "osmolyzed": "negative", "osmolyzing": "negative", "pair": "unknown",
    "paired": "unknown", "pairing": "unknown", "phosphorylate": "unknown", "phosphorylated": "unknown", "phosphorylating": "unknown",
    "phosphorylation": "unknown", "prevent": "negative", "prevented": "negative", "preventing": "negative", "prevention": "negative",
    "prohibit": "negative", "prohibited": "negative", "prohibiting": "negative", "prohibition": "negative", "produce": "positive",
    "produced": "positive", "producing": "positive", "production": "positive", "promote": "positive", "promoted": "positive",
    "promoting": "positive", "promotion": "positive", "react": "unknown", "reacted": "unknown", "reacting": "unknown",
    "reaction": "unknown", "reduce": "negative", "reduced": "negative", "reducing": "negative", "reduction": "negative",
    "regulate": "unknown", "regulated": "unknown", "regulating": "unknown", "regulation": "unknown", "relate": "unknown",
    "related": "unknown", "relating": "unknown", "relation": "unknown", "release": "unknown", "released": "unknown",
    "releasing": "unknown", "repress": "negative", "repressed": "negative", "repressing": "negative", "repression": "negative",
    "restrict": "negative", "restricted": "negative", "restricting": "negative", "restriction": "negative", "silence": "negative",
    "silenced": "negative", "silencing": "negative", "slice": "negative", "sliced": "negative", "slicing": "negative",
    "stimulate": "positive", "stimulated": "positive", "stimulating": "positive", "stimulation": "positive", "stop": "negative",
    "stopped": "negative", "stopping": "negative", "strap": "unknown", "strapped": "unknown", "strapping": "unknown",
    "subdue": "negative", "subdual": "negative", "subdued": "negative", "subduing": "negative", "supplement": "positive",
    "supplementation": "positive", "supplemented": "positive", "supplementing": "positive", "suppress": "negative", "suppressed": "negative",
    "suppressing": "negative", "suppression": "negative", "tether": "unknown", "tethered": "unknown", "tethering": "unknown",
    "trigger": "positive", "triggered": "positive", "triggering": "positive", "ubiquitination": "unknown", "ubiquitinylate": "unknown",
    "ubiquitinylated": "unknown", "ubiquitinylating": "unknown", "unite": "unknown", "united": "unknown", "uniting": "unknown",
    "union": "unknown", "weaken": "negative", "weakened": "negative", "weakening": "negative", "wrap": "unknown",
    "wrapped": "unknown", "wrapping": "unknown", "xylosylate": "unknown", "xylosylated": "unknown", "xylosylating": "unknown",
    "xylosylation": "unknown", "zip": "unknown", "zipped": "unknown", "zipping": "unknown",
}

def __label_relationship(BERTOutput: pd.DataFrame):
    '''
    Labels the relationships in the given DataFrame based on predefined term relationships.
    '''

    # Build data.frame of term relationships
    TermTable = pd.DataFrame(term_relationships.items(), columns = ["Term", "Relationship"])

    # Make a holder for all terms and descriptions
    all_terms = []
    all_descriptors = []

    for sentence in BERTOutput["Sentence"]:

        # Now IMPORTANT - limit the terms to the range between the two entities in the sentence
        sentence = sentence.split("@TERM$1")[1].split("@TERM$2")[0]
        
        # Extract the terms
        the_terms = [term for term in TermTable["Term"] if term in sentence]

        # If the length is more than 1, label each term
        if len(the_terms) > 1:
            
            # Extract desciptors
            the_descriptor = [TermTable[TermTable["Term"] == term]["Relationship"].values[0] for term in the_terms]

            # Collapse input
            the_terms = "; ".join(the_terms)
            the_descriptor = "; ".join(the_descriptor)

        else:
            the_terms = "No term detected"
            the_descriptor = "No term detected"

        # Add terms and descriptors
        all_terms.append(the_terms)
        all_descriptors.append(the_descriptor)

    # Add terms and descriptors 
    BERTOutput["Terms"] = all_terms
    BERTOutput["Descriptors"] = all_descriptors 

    # Return these descriptors as a new column in the BERT output, which can be used for further analysis.
    return BERTOutput