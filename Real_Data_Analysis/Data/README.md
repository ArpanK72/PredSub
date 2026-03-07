# Data

This folder contains the raw data files used in the real data analysis.
```
Data/
├── DBLP/
│   ├── coauth-DBLP-full-nverts.txt
│   ├── coauth-DBLP-full-simplices.txt
│   └── coauth-DBLP-full-times.txt
├── Cannes/
│   └── Cannes2013_multiplex_edges.txt
└── README.md
```

---

## DBLP Co-Authorship Network

### Source
Download from the Austin R. Benson hypergraph datasets repository:
[https://www.cs.cornell.edu/~arb/data/coauth-DBLP/](https://www.cs.cornell.edu/~arb/data/coauth-DBLP/)

### Description
The DBLP co-authorship dataset represents academic publications as hyperedges, where each hyperedge connects all co-authors of a publication. The dataset includes timestamps (years) for each publication. `DBLP_Preprocessing.R` converts the hypergraph into two co-authorship networks for Era 1 (2011--2014) and Era 2 (2015--2018), restricted to authors active in both eras with degree at least 5 in both networks.

### Citation
A. R. Benson, R. Abebe, M. T. Schaub, A. Jadbabaie, and J. Kleinberg. Simplicial closure and higher-order link prediction. *Proceedings of the National Academy of Sciences*, 2018.

---

## Cannes 2013 Twitter Multiplex Network

### Source
Download from the Manlio De Domenico multiplex network datasets repository:
[https://manliodedomenico.com/data.php](https://manliodedomenico.com/data.php)

### Description
The Cannes 2013 dataset captures Twitter interactions among users during the Cannes Film Festival in 2013. `Cannes_Preprocessing.R` extracts the Retweet and Mention layers, converts both to undirected networks on a common node set, and removes isolated nodes.

### Citation
E. Omodei, M. De Domenico, and A. Arenas. Characterizing interactions in online social networks during exceptional events. *Frontiers in Physics*, 3:59, 2015.
