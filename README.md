### Overview 
**Chemical fingerprints encode mother-offspring similarity, colony membership, relatedness and genetic quality in fur seals**  
Stoffel, M.A., Caspers, B.A., Forcada, J., Giannakara, A., Baier, M.C., Eberhart-Phillips,
  L.J. , Müller, C. & Hoffman, J.I.

**analysis_markdown.pdf** contains the code for our analyses. 
In addition, we provide our data files, which should be placed in a subfolder of the working directory called *files*.

* **genotypes.txt** contains raw microsatellite data with each locus in two columns and individuals in rows
* **scent_raw.csv** contains aligned GC data
* **factors.csv** contains data about colony membership (colony), mother-offspring pairs (family) and mothers and pups, respectively (age)
* **simper_mp_results.csv** contains average contributions of compounds in chemical fingerprints to mother-offspring similarity (output from Primer-E)
* **coordinates_beach1.csv** contains euclidian coordinates for seals on the special study beach
* **bootstrap_mums.csv** contains results from the BIO-ENV bootstrap procedure (see Rmarkdown file or Methods of the paper)
* **relatedness.csv** contains pairwise relatedness values for all individuals (Queller&Goodnight)