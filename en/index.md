---
layout: default-e
title: SRATS2010
---
# SRATS: Software Reliability Assessment Tool on Spreadsheet

## Overview

SRATS is a Microsoft Excel AddIn for estimating software reliability with non-homogeneous Poisson process (NHPP) based software reliability growth models (SRGMs). The features of SRATS are

- Excel interface is used for data arrangement,
- 11 types of NHPP-based SRGMs are implemented,
- 3 kinds of data form are available; and
- stable parameter estimation is available.

The our parameter procedures provide maximum likelihood (ML) estimates for any patterns of data, and are based on our research results.

## Installation

SRATS consists of two files:

- srat2010.xla,
- SRATS2010.dll,

which are zipped into a file. The installation is the same as ordinary Microsoft Excel AddIns. That is, after extracting two files from the zip file and copy two files to user's `AppData\Roaming\Microsoft\AddIns` folder. Next, start Excel and select Add-Ins of Excel options and check SRATS.

## Download and system requirements

The zip file for SRATS can be downloaded from <a href="./SRATS2010_en_20130610.zip">here (SRATS2010_en_20130610.zip)</a>

We recommend it runs on Microsoft Excel 2010 (32bit). If an error appears on the starting SRATS, please execute the following command on command prompt with Administrator:

```
regsvr32 C:\Windows\System32\MSCOMCTL.OCX # for 32bit OS
regsvr32 C:\Windows\SysWOW64\MSCOMCTL.OCX # for 64bit OS
```

## Quick start

SRATS treats three data types. Here we introduce grouped data which is most frequently used. The grouped data is given by the following two columns

|time interval|faults|
|:---:|:---:|
|3 | 5 |
|2 |  1|
|3.5|2|
|...|...|

The first and second columns indicate the time interval for observation and the number of discovered faults in the observation. Then the first row in the above example means 5 faults were discovered by 3 days testing.

After selecting the range representing grouped data (without headers), we start SRATS from AddIn menu. The 'range' on the main form of SRATS indicates the selected data. On the main form, push the `Estimate All` button for estimating model parameters for 11 models with the selected data. After finishing for parameter estimation, the list shows AIC (Akaike information criterion) and BIC (Bayesian information criterion) for all the model. Select the model with the smallest AIC or BIC on the list. Push `report...` and then the form is opened to draw and output the detailed results to new sheets. The click on OK in the form, SRATS draws the mean value function of the selected NHPP-based model with Excel functions.

Have fun!

## Contacts

If there is any trouble, contact to `okamu [at] rel.hiroshima-u.ac.jp`.
