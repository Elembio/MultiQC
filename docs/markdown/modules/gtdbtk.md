---
title: GTDB-Tk
displayed_sidebar: multiqcSidebar
description: >
  Assigns objective taxonomic classifications to bacterial and archaeal genomes.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/gtdbtk/gtdbtk.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Assigns objective taxonomic classifications to bacterial and archaeal genomes.

[https://ecogenomics.github.io/GTDBTk/index.html](https://ecogenomics.github.io/GTDBTk/index.html)
:::

The module parses `summary.tsv` outputs from GTDB-Tk's `classify.py` and `classify_wf.py`.

The module only works for version >= 2.4.0 because column names changed.

`classify.py` and `classify_wf.py` are used to determine the taxonomic classification of input genomes.

### File search patterns

```yaml
gtdbtk:
  contents: "user_genome\tclassification\tclosest_genome_reference\tclosest_genome_reference_radius\t\
    closest_genome_taxonomy\tclosest_genome_ani"
  num_lines: 10
```
