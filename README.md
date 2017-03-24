# SNPversity #

Web-based tool to visualize single nucleotide polymorphisms (SNP) diversity.
Can be accessed [here](http://www.maizegdb.org/snpversity).

## Structure ##

#### [/time_estimate](/time_estimate)
This folder contains the [python script](/time_estimate/fetch_time.py) used to predict query execution time. It accepts dataset, # stocks, range of positions (in bp) as an input, and returns the estimated time in seconds. For example, in order to get estimated query processing time for AllZeaGBSv2.7 dataset with 100 stocks across 18121 bp's:

```
python2.7 fetch_time.py AllZeaGBSv27public20140528 100 128121
```

#### [/tassel](/tassel)
Description here

#### Folder 3
Description here


## Environment Requirements ##

#### Java Version 1.8.0_77
* Dependency 1
* Dependency 2

#### PHP Version 5.3.3
* No additional libraries required.

#### Python Version 2.7.6
* See [requirements.txt](requirements.txt).

#### Database (PostgreSQL)