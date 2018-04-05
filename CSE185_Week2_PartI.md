# Week 2: Why did I get the flu? (part I)
Skills covered: review of tools from week 1, deep sequence data, awk one-liners, git

This year, you prudently got the flu vaccine. So when your roommate, who forgot to get vaccinated, came down with the flu, you weren’t worried. But somehow, several days later, you started feeling feverish, weak, and all-around-awful. You knew it was the flu, but how could this be?!

Suspecting that this season’s vaccine wasn’t a match for the flu virus that infected your roommate, you had some friends in the med school run a hemagglutination inhibition (HI) assay on virus samples from your roommate.

The results showed that your roommate’s virus closely matched the HI profile for an H3N2 strain called A/Hong Kong/4801/2014 (H3N2). Was that one of the flu strain’s covered by this season’s (2017/2018) vaccine? Find out what strains were in this season’s vaccine, then answer the IClicker question.

You’ve heard of viral quasispecies, and suspect that maybe a small portion of the virus population mutated and evolved while replicating inside your roommate’s cells, which could explain how it was able to infect you. To find out, you have your friends set up a targeted deep sequencing experiment to analyze the HA genes in your roommate’s viral sample. They set up an Illumina single-end sequencing run. When they send you the results, you start analyzing your roommate’s sequence right away.

## 1. Git from the command line
TODO setting up lab notebook and report from the command line

## 2. Inspect the data from your roommate

This sequencing data can is in the `public/week2` directory. Your roommate’s data is labeled
SRR1705889. Record how many reads there are in this file, then look at the first 20 lines with the
head command and answer the IClicker question.

Reminder - to get from your home directory to the parent cs185s class directory:
```shell
cd ..
```

Since there seem to be reads of multiple lengths, you suspect the data may have been processed
somehow. You can use a bash one-liner with `awk`, `sort`, and `uniq` to find out. Below is a reminder for
the general syntax of awk.

```shell
awk '/search_pattern/ {actiontotakeonmatches; otheractions;}' file_to_awk
```

For the actual command (below), use pipes (`|`) to send the output of one command as input to the next. 

```shell
cat SRR1705889.fastq | awk 'NR%4==0 {print length}' | sort -n | uniq -c
```

Let's walk through the commands we piped together:
* `cat`: This command, which stands for "concatenate", just outputs the contents of a file. Here we are opening the SRR file with `cat` and piping it into `awk`, so we don’t have to specify the file to `awk` after the brackets. 
* `awk`: `NR%4`, means NR ("number of records", or line number) modulo 4. The command returns every 4th line (because 4 goes into 4, 8, 12, 16, etc. perfectly with no remainder). This gives us a quick way to scan through each read (remember each read has 4 lines). The 4th line is the quality string, with a value for each base, so its length is the same as the length of the read.
* `sort` sorts the output lines. The `-n` flag tells it it's dealing with numbers so sort numerically. The next command requires the input to be sorted first.
* `uniq` returns only unique lengths (many reads are probably the same length), and the `-c` flag gives a count of how many there are of each.

**Record** the maximum read length, then answer the IClicker again. It looks like the reads come in various lengths, which suggests that someone has already processed the data. You check
with your friends in the med school, and it turns out that yes, they trimmed the low quality bases from
the ends of the reads for you.

## 3. Align your roomates data to the reference sequence

The reference sequence for the influenza hemagglutinin gene is not in the public folder, so you will have to
download it from NCBI. This can be done with a command from the EntrezDirect utility, which is
installed for this course. 

`cd` into your working directory, make a folder for `week2`, and use the EntrezDirect command `efetch` to download the reference sequence, which has the NCBI id number: KF848938.1. An example of the `efetch` usage is below, you can also type `efetch -help` for more information. You must specify which NCBI database to use, the id or ascension number of the sequence, and the format you would like. Redirect the output into a file with the “.fasta” extension:

```shell
efetch -db nucleotide -id KF848938.1 -format fasta > KF848938.1.fasta
```

Index the reference file with bwa:
```shell
bwa index KF848938.1.fasta
```

Align your roommate’s viral data to the reference sequence and make an mpileup. Since we learned
how this works last week, you can do most of it in one line today with pipes. This will also avoid
making all those intermediate files. 

```shell
bwa mem KF848938.1.fasta /pathto/SRR1705889.fastq | \
    samtools view -S -b | \
    samtools sort > roommate.bam
```

Run samtools view with the `-f 4` command to extract all unmapped reads, then count them. Use the number of unmapped reads and the total reads you calculated from the unprocessed `SRR1705889.fastq` file to **calculate the number of reads that mapped.** 

```shell
samtools view -f4 roommate.bam | wc -l
```

For a nice explanation of SAM flags, see this tool: broadinstitute.github.io/picard/explain-flags.html. 

Finally, **index** the bam file:
```shell
samtools index roommate.bam
```

## 5. Look for common variants with VarScan

Make an mpileup of the bam alignment file. To save computing power, `samtools` default behavior
stops piling up the base calls at each position when it gets to 8000 calls. Since our variants may be
quite rare, set that depth limit to something we know is higher than our coverage with the `-d` flag.

```shell
samtools mpileup \
    -d 1000000 \
    -f KF848938.1.fasta \
    roommate.bam > roommate.mpileup
```

Run VarScan on the mpileup. You may have to `cp` the VarScan.jar file from `week1` to `week2`. First, look for positions where most of the viruses infecting your roommate differ from the reference. Maybe there is a common mutation that wouldn’t have shown up in the HI test. Use a high minimum variant frequency cut-off to find only those mutants present in most (95% or more = 0.95) of the viral DNA molecules. 

```shell
java -jar VarScan.jar mpileup2snp \
    roommate.mpileup \
    --min-var-freq 0.95 \
    --variants --output-vcf 1 > roommate.vcf
```

How many variants are reported back? To clean up the data and allow you to easily copy and paste,
**edit the awk script below so it only output the fields you are interested in** (position, reference base,
alternate base, in that order). Run the command, then **record the reference base, position, and mutant
base in your lab notebook.**

```shell
cat roommate.vcf | awk 'NR>24 {print $1, $2}' # Note you'll need to edit this
```

What do these mutations do? Could they be what allowed your roommate’s virus to escape the antibodies in your body from the flu vaccine? Since we are only looking at a single gene (and since the VEP doesn’t have a function for viruses) we will do this manually. In a browser, open the free, online sequence editor called WebDSV at
http://www.molbiotools.com/WebDSV/index.html

Go to NCBI (in a browser) and search for the reference sequence KF848938.1. Click on the FASTA link in the first result, then copy and paste the DNA sequence (not the header) where it says ‘paste sequence here’ on WebDSV. Click the yellow "process" button, then click "select all" and "translate". Close the pop up window, and when you get back to WebDSV, you should see the single-letter amino acid code above your sequence in the bottom window.

In the right side of the app, click the button with "AA" on it to open a codon table. For each position in your list, find the codon it is a part of. (You can hover over the bases to show the position, and you can hover over the amino acids to show the codon). Use the codon table to translate each mutation, and **record the original codon, mutated codon, original amino acid, its position in the protein, and the mutated amino acid. Then, record
whether the change is synonymous or non-synonymous.**

Example:
A72G ACA>ACG Thr24Thr synonymous

Repeat with the rest of the variants, record the results, and answer the IClicker question.

## 4. Look for rare variants with VarScan

Now try looking for rare variants. Set the min var freq to 0.001 (0.1%) and run the scan again, on the
same mpileup file. 

```shell
java -jar VarScan.jar mpileup2snp \
    roommate.mpileup \
    --min-var-freq 0.001 \
    --variants --output-vcf 1 > roommate_rare.vcf
```

How many variants are reported back now, and how abundant are they? Modify the `awk` script again
to extract the position, reference base, and frequency. The frequency is part of field `$10`, so run the
command below.

```shell
cat roommate_rare.vcf | awk 'NR>24 {print $2, $4, $5, $10}'
```

This DOES contain the frequency, but it’s buried in that last really long field. Fortunately, awk can use
more than just spaces as field separators. Pipe this output into the awk command below to process it
further. The `-F "[xy]"` option, splits it into fields separated by EITHER x OR y. So now it sees our 3
original fields (separated by " "), plus all the additional fields separated by ":" . Edit the command
below to output the field with the frequency (the \%), and then run the whole pipeline starting
from `cat`.

```shell
awk -F '[ :]' '{print $1, $2, $3, $? }'
```

Copy the results into your notebook. 

You take your data back to your friends at the med school and excitedly show them all of the rare
variants you found. They don’t seem nearly as excited as you, and when you ask why, they point out
that at these frequencies its very difficult to tell the difference between real, rare mutants in the viral
population, and errors introduced in the sequencing and amplification process. They suggest a control reaction, where they will sequence an isogenic viral sample (all virus particles genetically identical) derived from a virus clone that matches the reference sequence. By looking at the errors in this reference, and comparing them to the mutations in your roommate’s sample, you hope you’ll be able to figure out which variants are real. Next time, we’ll analyze this control sequencing data. 

**Acknowledgements**: Adapted from a lab originally written by Dr. Katie Petrie
