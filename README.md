# Drosophila genomics analysis

## To-do
- [ ] Move unorganized scripts
- [ ] Clean code on unorganized scripts
- [ ] There are universal functions inside 03_cmh_tests.R that have to be moved somewhere else
- [ ] Unstage big files
- [ ] MAF excluding B populations
- [ ] Allele trajectory means instead of one per replicate
- [ ] Organize my code
- [ ] Vectorize function calling
- [ ] Scaled Classic CMH permutation test

## Issues
 - Biggest issue I'm facing - my CMH test p-values are too small! When I check how other papers in the field are dealing with FDR, lots of them apply a correction such as q-value or Bonferroni or BH, and then assume everything < 0.05 as true positives. However, if I do that, half of my SNPs are significant, or in other words, the whole Drosophila genome becomes associated with the adaptation.

I can go to [this page](https://scholar.google.com/scholar?hl=en&as_sdt=5%2C38&sciodt=0%2C38&cites=4613566741279609699&scipsc=1&q=drosophila&btnG=) and see all the papers that cited the ACER (i.e. all the papers that cited the adapted.cmh.test).
I filter by Drosophila only and I got a few matches.
[This paper is interesting](https://www.biorxiv.org/content/10.1101/2024.02.08.579525v1.full.pdf), as it also appears to have this problem, and they simply apply a FDR threshold of 0.01 (it's not clear if its p-value = 0.01 or if they are assuming their top 0.01 SNPs as meaningful).
[This](https://academic.oup.com/genetics/article/224/3/iyad050/7085646) is also a good paper, but their p-values are higher.
[This](https://link.springer.com/article/10.1186/s13059-021-02425-9) is a great paper with higher p-values.
[This](https://www.nature.com/articles/s41467-022-31622-8#data-availability) is an excellent paper, that uses copepods and also have higher p-values. Their code is very messy, though.

# 2025-02-24
I set my significance level to be 100. That is an extremely high number and I think we will get in trouble with reviewers because of this. However, right now, I can't think of a better solution given the infinitesimal p-values we have.

My code is mostly clean, but not quite.