# Deep Phenotyping
This study was a collaboration I worked on between the [Buckner lab](https://cnl.rc.fas.harvard.edu/) at Harvard's Center for Brain Sciences and the [Onnela lab](https://www.hsph.harvard.edu/onnela-lab/) at the Harvard T.H. Chan School of Public Health. Large amounts of data were collected from a small number of individuals to develop methods of deeply characterizing behavior at the individual level. This included activity and sleep estimates via wrist-worn monitors, regular surveys sent via smartphone, and passively collected smartphone data such phone movement or whether the phone screen was on. From this information, we were able to deeply characterize individuals’ behaviors and changes in those behaviors over time.

### social_pca.R
Using survey response (e.g. how lonely did you feel today?) along with text and call data, I conducted a factor analysis to examine latent constructs underlying social behavior.

### sleep_prediction_imp_pmm.R
Using only data passively collected from smartphones, I trained a logistic regression to predict with 95% accuracy whether a person was sleeping (assessed by wrist-worn monitors) at any given minute.

