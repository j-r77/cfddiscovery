# Discovering Conditional Functional Dependencies

> Based on the conference paper "Revisiting Conditional Functional Dependency Discovery: Splitting the “C” from the “FD”.", by J. Rammelaere and F. Geerts, to be published in the Joint European Conference on Machine Learning and Knowledge Discovery in Databases (ECML PKDD 2018). Code by Joeri Rammelaere.

## Quick Start Guide
The code compiles with cmake. 
The algorithm takes four mandatory arguments and one optional as input:
1. A dataset in csv format
2. The minimum support threshold
3. The confidence leeway threshold
4. The maximum antecedent size of the discovered CFDs
5. [Optional] The algorithm/implementation to be used (See "Algorithmic Choices")

## Available Datasets
The Data/ folder contains the 3 datasets used in the experimental section, as well as the Contraceptive dataset, taken from the UCI Machine Learning Repository (http://archive.ics.uci.edu/ml/). 

Statistics of the datasets:

Dataset | Nr. Tuples | Nr. Items | Nr. Attributes
:--- | :--- | :--- | :---
Adult | 48842 | 202 | 11
Mushroom | 8124 | 119 | 15
Nursery | 12960 | 32 | 9

## Problem Statement
Given an instance D of a schema R, support threshold delta, confidence leeway threshold epsilon, and maximum antecedent size alpha, the approximate CFD discovery problem is to find all CFDs phi: (X -> Y, tp) over R with:
* support(phi,D) >= delta
* confidence(phi,D) >= 1-epsilon
* |X| <= alpha.

## Algorithmic Choices
As described in the paper, the code supports 10 distinct algorithms, belonging to three different general methodologies: integrated, itemset-first, and fd-first. The integrated method can be run depth-first or breadth-first (corresponding to normal CTane). The itemset-first and fd-first methods consist of 2 steps, each of which can also be run either depth-first or breadth-first. 
This gives the following options:
* Integrated-BFS (CTane)
* Integrated-DFS
* Itemset-First-BFS-bfs
* Itemset-First-BFS-dfs
* Itemset-First-DFS-bfs
* Itemset-First-DFS-dfs
* FD-First-BFS-bfs
* FD-First-BFS-dfs
* FD-First-DFS-bfs
* FD-First-DFS-dfs

By default, the option FD-First-DFS-dfs is chosen, as it is typically the fastest.

## Paper Abstract
Many techniques for cleaning dirty data are based on enforcing
some set of integrity constraints. Conditional functional dependencies
(CFDs) are a combination of traditional Functional dependencies (FDs)
and association rules, and are widely used as a constraint formalism for
data cleaning. However, the discovery of such CFDs has received limited
attention. In this paper, we regard CFDs as an extension of association
rules, and present three general methodologies for (approximate) CFD
discovery, each using a different way of combining pattern mining for discovering
the conditions (the “C” in CFD) with FD discovery. We discuss
how existing algorithms fit into these three methodologies, and introduce
new techniques to improve the discovery process. We show that the right
choice of methodology improves performance over the traditional CFD
discovery method CTane.
