# Week 2: Why did I get the flu? (part II)
Skills covered: estimating coverage, understanding sequencing errors

Realizing that you can’t tell the difference between sequencing errors and real mutations in the HA
gene from rare members of the viral quasispecies, you turn to your med school friends for help. They
suggest using a control reaction to better understand the sequencing errors.

After they sequenced your roommate’s sample, which you think contains a mixture of virus
genotypes, they next sequenced three control samples, consisting of isogenic (100% pure) samples
of the standard (reference) H3N2 influenza virus.

Any “mutations” you detect in the control samples, which don’t contain any true genetic variants, must
be due to sequencing errors. You can use the frequency of the sequencing errors from the control to
help figure out what’s an error and what’s a true variant in the data from your roommate.

*Note: Today, you will start seeing fewer commands pasted into the tutorial, especially for things we have
already done. You are more than welcome to go back to previous tutorial or to your notebook entries,
and you can also usually type a command (by itself or with `-help` or `-h`) to bring up the usage guide.
You can look up documentation for particular bits of software online, and you can definitely ask
the instructor, TA, or your neighbors for help if you get stuck or aren’t certain about a command you’ve
typed. Pay attention to the text that pops up when commands are running, as it can often give you a
clue about potential problems.*

## 6. Inspect and align the control sample sequencing data

Fastq data for the three controls (from sequencing of isogenic reference samples) are in the
`public/week2` folder. The three controls are called `SRR1705858.fastq`, `SRR1705859.fastq`, and
`SRR1705860.fastq`. Calculate how many reads are in each file, and **record it in your notebook.**

**Calculate a rough estimate of the coverage** in your samples. Since all 3 controls and your roommates
sample have about the same number of reads, you can use any one of the read counts for this. For
read length, use 151, even though we know some of the reads are shorter. Remember, coverage is
the number of reads corresponding to each position in the reference, so you also need to know how
many base pairs long the reference sequence is. Include the result and how you calculated it in the
lab report. 

`cd` into your working directory, and use the `bwa mem` pipeline from last time to align each control
fastq file to the reference (`KF848938.1.fasta`), convert it to a bam file, and sort it. You need to
make a separate bam alignment for each control sample, but you do not have to index the reference
again. Be sure to give each bam file a unique name. Use samtools to index each of the new control alignment (bam) files. 

**Automating the analysis (optional, but highly recommended)** To make your life easier, you may want to take advantage of bash variables and for loops so you don't have to copy the same code over and over again. Below is an example of how to set and access a variable in bash:

```
myVariable=200
echo $myVariable
```

Note, you must not include any extra spaces around the "=" and you need to use a $ symbol to access the variable. `echo` simply prints results to the terminal's standard output.

If you'd like to run the same command multiple times but slightly altered, you can use a for loop. e.g.:
```
for x in sample1 sample2 sample3
do
  echo "let's run some commands on $x"
  echo "bwa mem /path/to/index /path/to/$x.fastq
done
```
(Of course, you will need to remove the `echo` and edit the paths in the `bwa` command above).

<blockquote>
**UNIX TIP**: If you want to run a multiple line command, you can put all the lines separated by semicolons and run on one command. e.g. "for x in sample1 sample2; do echo $x; done". Alternatively, you can enter the commands line by line or copy and paste the entire multi-line command from a text file such as your lab notebook.
</blockquote>


## 7. Use VarScan to look for rare variants in the reference files.

Run VarScan, with a minimum variant frequency of 0.001 (0.1%) on each of the reference alignments. Be sure to tell VarScan to only output variants, and to format the output in the vcf format.

Use awk to parse each vcf file so that you get three lists containing position, reference base, alternative base, and frequency. Copy those lists into a spreadsheet. You should upload this spreadsheet and cite as supplementary data in your lab report. 

## 8. Compare the control results to your roommate results

Examine the data in the spreadsheet. Are there any positions in your roommate’s sample (aside from the common ones from last week) that stick out as possibly being more than just sequencing error? Calculate the average and standard deviation of the frequencies reported within each reference (you should have 3 averages and 3 standard deviations, and put them in the same units). 

Answer the IClicker about the error rate of sequencing.

Did VarScan report rare mutations in your roommate’s file with frequencies that are more than 3 standard deviations away from the averages in the reference files? If there are, use WebDSV to identify the original amino acid at each position, position number in the protein sequence, and the amino acid resulting from the mutation. Record these in your notebook. 

Discussion question: Are there any positions reported by VarScan in all 3 of the reference
sequences? You could, in principle, also calculate the average and standard deviation between the 3
reference replicates for one position at a time. Which kind of average and standard deviation do you
think is better for error correction?

## 9. Epitope mapping

Use the epitope locations listed in Munoz et al (listed under reading for the lab) to determine if any of the high confidence (> 3 std deviations away from reference error rate) mutations from your roommate’s flu infection are located in an epitope region of hemeagglutinin. (Epitopes are the parts of the protein structure recognized by antibodies). If so, list which epitope regions are mutated.

## 10. Lab Report

For the abstract, remember that your goal is to figure out how you got the flu even though the HI test
said your roommate’s flu strain was covered by the vaccine.

For the introduction, briefly cover how the flu vaccine works, and the idea of antigenic drift and viral
quasispecies. Also provide some background on targeted deep sequencing for studying mixed
populations, and the sources of error in next generation sequencing that this lab project tried to
correct.

In your methods section, briefly introduce the sequencing data (what are the samples, how many
cycles was the sequencing run, was the data pre-processed?). Describe how you used VarScan to
examine both your roommate’s sample and the reference samples. Be sure to describe the
parameters you used if they were not the default options.

In your results section, include the number of reads you started with and the number of reads that
mapped, for each of the 4 data sets. Report the average and standard deviation of the frequencies
from each reference sample. For your roommates data, make tables or a list containing the common
mutations, and any rare mutations you think are truly represented in the viral population. List these
mutations both in terms of the DNA base, and the protein amino acid changes, and state whether
they are non-synonymous or synonymous. You do not need to list all of the mutations reported by
VarScan in your lab report, but include them in your notebook.

State whether any of the mutations affect the epitope regions of the hemeagglutinin protein and if
so, report which epitopes they affect.

In the discussion, be sure to explain how you decided which mutations were most likely to be real.
Given your results, explain how you think you were able to get the flu from your roommate, even
though you had received the flu vaccine. Also in the discussion section, propose three additional
ways to control for error in deep sequencing experiments like this, and explain why error control is
important for accurately identifying and quantitating rare variants. Our approach was pretty quick and
simple; there are many more sophisticated methods out there, some of which we talked about in
class. You can suggest laboratory steps to minimize errors in the first place, bioinformatics steps you
could implement on our data, or existing software. For each suggestion, include a sentence or two
explaining how it would reduce error. 

**Optional Extra-Credit Challenge Question** For up to 1 bonus point on this lab report. How would
you calculate the ACTUAL average coverage for one of our data sets, only for mapped reads, and
taking into consideration the fact that the reads are not all the same length? If you use a script or
software that someone else wrote, please explain how it works, and how you would call it at the
command line. Include your approach, and your answer, if you found one, at the end of your lab
report, after the discussion. 
