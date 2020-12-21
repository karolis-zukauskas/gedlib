# build
  - python3 install.py --debug --lib gxl --gurobi /mnt/d/Masters/gurobi910 --test all
  - python3 install.py --debug --lib gxl --gurobi /mnt/d/Masters/gurobi910 --test kazu

# Apatinis rezis
Best:
  - [LP] ADJ-IP, F2
  - [Mixed] BRANCH-TIGHT

Fastest:
  - [LSAP] NODE, BRANCH-FAST, BRANCH-CONST

Optimizuoti **ADJ-IP (BLP_NO_EDGE_LABELS)** ??


# Virsutinis rezis
Best:
  - [LS] IPFP

Fastest:
  - [LSAP] NODE, BRANCH-FAST, BRANCH-CONST, RING, STAR

Optimizuoti **ADJ-IP** ??