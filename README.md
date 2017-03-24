# SNPversity #

Web-based tool to visualize single nucleotide polymorphisms (SNP) diversity.
Can be accessed [here](http://www.maizegdb.org/snpversity).

## Structure ##

#### [/time_estimate](/time_estimate)
This folder contains the [python script](/time_estimate/fetch_time.py) used to predict query execution time. It accepts dataset, # stocks, range of positions (in bp) as an input, and returns the estimated time in seconds. For example, in order to get estimated query processing time for AllZeaGBSv2.7 dataset with 100 stocks across 18121 bp's:

```
$ python2.7 fetch_time.py AllZeaGBSv27public20140528 100 128121
```

#### [/tassel](/tassel)
The [Tassel](http://www.maizegenetics.net/tassel) part of the project used for processing queries lives here. The Tassel wrapper source code is contained in file [Tassel_gt_server.java](/tassel/tassel-wrapper/src/tassel_gt_server/Tassel_gt_server.java). Other interesting folders are:

* [tassel-gt-server-david-SPLIT.jar](tassel/tassel-gt-server-david-SPLIT.jar) is the exported JAR file which gets called by the server-side PHP scripts.

#### Folder 3
Description here


## Environment Requirements ##

#### Java Version 1.8.0_77
* The Tassel files (version 5.2.17) are included in the repository for your convenience. All other required libraries need to be added.
* For a list of these libraries and where to put them, consult [manifest.txt](/tassel/tassel-wrapper/manifest.txt)

#### PHP Version 5.3.3
* No additional libraries required.

#### Python Version 2.7.6
* See [requirements.txt](requirements.txt).

#### Database (PostgreSQL)