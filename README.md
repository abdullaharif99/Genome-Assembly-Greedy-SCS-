# Genome-Assembly-GreedySCS
A python script for genome assembly using the greedy shortest common superstring algorithm

# Files

The input is provided in the form of mystery_virus.fq file.
All the reads in mystery_virus.fq contain 100 base pairs each.

The script outputs two file: output_stats.txt which contains some basic information regarding the assembled sequence and mystery_virus_genome.txt which contains the actual assembled genome.

# Usage
To use this script, just change the name of the input FastQ file in Line 48 to your required filename.

```python
reads = readFastq('your_filenmae.fq')
```

You can also change the output name in Line 51.

```python
with open('your_output_filename.txt', 'w') as f:
```
You can also adjust the minimum overlap in Line 49 to experiment with the results

```python
assembled_ss = greedy_scs(reads, k=your_min_length)
```

