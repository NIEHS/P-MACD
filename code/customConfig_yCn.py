motifs2Find = ("A", "T", "G", "C", "Cg", "cG", "[atgc]Tt", "aA[atgc]", "[ct]C[atgc]", "[atgc]G[ga]", "[atgc]Cg", "cG[atgc]")
findTitles = ("A", "T", "G", "C", "Cg", "cG", "nTt", "aAn", "yCn", "nGr", "nCg", "cGn")
motifs2Count = ("a", "t", "g", "c", "cg", "[atgc]tt", "aa[atgc]", "[ct]c[atgc]", "[atgc]g[ga]", "[atgc]cg", "cg[atgc]")
countTitles = ("a", "t", "g", "c", "cg", "ntt", "aan", "ycn", "ngr", "ncg", "cgn",)
countTitles = tuple([x+"_counts" for x in countTitles])
