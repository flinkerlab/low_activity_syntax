The code for this "Low Activity Syntax" project is broken down into 4 directories: 
 - "functions": some custom functions used throughout the scripts;
 - "preprocessing": scripts for a preprocessing a sample subject (these were largely the same across subjects, except for hard-coded trial and electrode exclusions);
 - "basic analyses": temporal warping, and sentence-list comparisons;
 - "RSA": scripts for all analyses involving RSIs

A brief glossary, translating from manuscript terminology to variable names in the code:
 - "RSI" (manuscript) = "zs.adjusted" (code)
 - "Syntax" (i.e., RSI type) = "diff.voice"
 - "Semantics" (i.e., RSI type) = "diff.event.semantics"
 - "(Sub)lexical" (i.e., RSI type) = "diff.1st.word"
 - "Sentence Production" (block/trials/task) = "sp"
 - "List Production" (block/trials/task) = "lp"
 - "Picture Naming" (block/trials/task) = "pn"

