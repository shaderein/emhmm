========================================================================

  Eye-Movement analysis with Hidden Markov Models (HMMs)
  emhmm-toolbox

  Antoni B. Chan, City University of Hong Kong
  Janet H. Hsiao, University of Hong Kong
  Tim Chuk, University of Hong Kong
  Cynthia Y.H. Chan, University of Hong Kong
  Lan Hui, City University of Hong Kong

  Copyright (c) 2021, City University of Hong Kong & University of Hong Kong

========================================================================

--- DESCRIPTION ---
This is a MATLAB toolbox for analyzing eye movement data with hidden Markov Models (HMMs). 
The major functions of the toolbox are:
  1) Estimating HMMs for an individual's eye-gaze data.
  2) Clustering individuals' HMMs to find common strategies.
  3) Visualizing individual and group HMMs.
  4) Statistical tests to see if two HMMs are different.
  5) Using representative Holistic and Analytic models to compute HA scale value.
  6) Computing entropy of HMMs
  7) Co-clustering individuals' HMMs (from different stimuli) to find common strategies
  8) Visualizing co-clustering group HMMs.
  
--- BRIEF INSTRUCTIONS ---
1) In MATLAB, run "setup" to setup the toolbox.
2) Note: if you are using Mac, you may need to disable Gatekeeper security for
   the MEX files in the toolbox. See Sec 3.1 in the documentation for more info.
3) In demo folder, run "demo_faces" for an example.
4) Also see "demo_faces_jov_clustering", "demo_faces_jov_compare", "demo_faces_duration", 
   and "demo_conversion".
5) There are also demos in the "models" folder: "demo_PBR_model".
6) In "demo" folder run "demo_vb_faces" and "demo_vb_faces_duration" for demos using 
   VBHEM for clustering, which automatically determines the number of clusters and states.
7) In "demo_cocluster" folder, run "demo_cocluster" for an example of using co-clustering. 
   The script "brm_cocluster_analysis" will analyze the co-clustering results from the BRM
   journal paper.

More documentation and descriptions are in the "docs" folder. In particular, 
"docs/emhmm-tutorial.pdf" contains slides introducing how to use the toolbox, while
"docs/emhmm-documentation.pdf" contains more detailed information.

--- REFERENCES ---
If you use this toolbox, please cite the following papers:

For learning HMMs for eye gaze data:
  Tim Chuk, Antoni B. Chan, and Janet H. Hsiao.
  "Understanding eye movements in face recognition using hidden Markov models."
  Journal of Vision, 14(11):8, Sep 2014.

For clustering HMMs with the VHEM algorithm:  
  Emanuele Coviello, Antoni B. Chan, and Gert R.G. Lanckriet.
  "Clustering hidden Markov models with variational HEM".
  Journal of Machine Learning Research (JMLR), 15(2):697-747, Feb 2014.

If you use the representative Holistic/Analytic models and HA Scale: 
  Cynthia Y.H. Chan, Antoni B. Chan, Tatia M.C. Lee, and Janet H. Hsiao.
  "Eye Movement Patterns in Face Recognition are Associated with Cognitive 
   Decline in Older Adults".
   Psychonomic Bulletin & Review, 25(6):2200-2207, Dec 2018.

For clustering HMMs with VBHEM (selecting the number of clusters automatically):
  Hui Lan, Ziquan Liu, Janet H. Hsiao, Dan Yu, and Antoni B. Chan.
  "Clustering Hidden Markov Models With Variational Bayesian Hierarchical EM."
  IEEE Trans. on Neural Networks and Learning Systems, To appear 2021.

For co-clustering HMMs with VHEM or VBHEM:
  Janet H. Hsiao, Hui Lan, Yueyuan Zheng, and Antoni B. Chan.
  "Eye Movement analysis with Hidden Markov Models (EMHMM) with co-clustering."
  Behavior Research Methods, April 2021.

--- CHANGE LOG ---
see CHANGELOG.txt for version information.

--- CONTACT INFO ---
Please send comments, bug reports, feature requests to Antoni Chan (abchan at cityu dot edu . hk).

--- ACKNOWLEDGEMENTS ---

This research was supported by the Research Grant Council of Hong Kong SAR: 
#17609117, #17402814 and HKU 745210H for J.H. Hsiao; CityU 110513 and 
G-CityU109/14 for A.B. Chan. We also thank the HKU Seed Funding Programme 
for Basic Research (Project numbers 201311159131 and 201811159165) to J.H. 
Hsiao.  We also thank the Strategic Research Grant from City University of 
Hong Kong (Project No. 7005218) to A.B. Chan.
