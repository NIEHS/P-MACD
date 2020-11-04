# source("/home/klimczakl/projects/yeast/context/findMotif.py")

motifs2Find = ("A", "T", "G", "C", "Cg", "cG", "tC[at]", "[at]Ga", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tC[atc]", "[atg]Ga", "cC", "Gg", "[at][ag]C", "G[ct][at]", "Cc", "gG", "[at]A", "T[at]")
findTitles = ("A", "T", "G", "C", "Cg", "cG", "tCw", "wGa", "tCa", "tGa", "tCt", "aGa", "tC", "Ga", "tCh", "dGa", "cC", "Gg", "wrC", "Gyw", "Cc", "gG", "wA", "Tw")

# for countMotifs.R
motifs2Count = ("a", "t", "g", "c", "cg", "tc[at]", "[at]ga", "tca", "tga", "tct", "aga", "tc", "ga", "tc[atc]", "[atg]ga", "cc", "gg", "[at][ag]c", "g[ct][at]", "cc", "gg", "[at]a", "t[at]")
countTitles = ("a", "t", "g", "c", "cg", "tcw", "wga", "tca", "tga", "tct", "aga", "tc", "ga", "tch", "dga", "cc", "gg", "wrc", "gyw", "cc", "gg", "wa", "tw")

#countTitles = map(lambda x: x+"_counts", countTitles)
countTitles = tuple([x+"_counts" for x in countTitles])

apobecTitles = ("tC_mutation", "tC_mutation_to_G", "tC_mutation_to_T", "APOBEC_mutation", "APOBEC_mutation_to_G", "APOBEC_mutation_to_T")

