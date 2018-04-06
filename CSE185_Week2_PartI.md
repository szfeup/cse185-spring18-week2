# Week 2: Why did I get the flu? (part I)
Skills covered: review of tools from week 1, deep sequence data, awk one-liners, git, (optional) emacs

This year, you prudently got the flu vaccine. So when your roommate, who forgot to get vaccinated, came down with the flu, you weren’t worried. But somehow, several days later, you started feeling feverish, weak, and all-around-awful. You knew it was the flu, but how could this be?!

Suspecting that this season’s vaccine wasn’t a match for the flu virus that infected your roommate, you had some friends in the med school run a hemagglutination inhibition (HI) assay on virus samples from your roommate.

The results showed that your roommate’s virus closely matched the HI profile for an H3N2 strain called A/Hong Kong/4801/2014 (H3N2). Was that one of the flu strain’s covered by this season’s (2017/2018) vaccine? Find out what strains were in this season’s vaccine, then answer the IClicker question.

You’ve heard of viral quasispecies, and suspect that maybe a small portion of the virus population mutated and evolved while replicating inside your roommate’s cells, which could explain how it was able to infect you. To find out, you have your friends set up a targeted deep sequencing experiment to analyze the HA genes in your roommate’s viral sample. They set up an Illumina single-end sequencing run. When they send you the results, you start analyzing your roommate’s sequence right away.

## 1. Setting up
First of all, since you only have a limited (~10GB) disk space quota for the course, let's clean out our directories from last week. You won't be needing any of the data from Lab 1 (and even if you do, your lab notebook should have everything you need to reproduce your results!). So go ahead and use `rm` to remove data files. You can remove files one-by-one (`rm filename`), but that will take a while. Typing `rm week1/*` should remove all the files in your week 1 folder. Please be careful with this command! You should always triple check before running an `rm` command to make sure you don't irreversibly delete something you will need later.

Now, we'll first get set up using Git from the command line. Last week, we edited all the files from the web browser. But we can also get all these files to edit through the terminal. 

Git is a really useful framework for keeping track of changes to files both for yourself or when working with a team of people. It maintains files in a "repository". You can update the repository through "commits". Each time you commit to the repository, Git will keep track of exactly what changes were made along with a short message describing the purpose of the changes. It can get infinitely complicated but today we'll go through some basic functionality.

To get started, we will make a "clone" (copy) of your assignment repository in your home directory in `ieng6.ucsd.edu`. Copy the base URL for the your repository (e.g., https://github.com/cse185-sp18/cse185-week2-<username>/). Log into the cluster and navigate to your home directory for the course using `cd`. Clone the repository using the following command:
```
git clone https://github.com/cse185-sp18/cse185-week2-<username/ week2
```

This will make a copy of the entire repository in a folder named `week2`. We can now edit and add files to the repository directly from the command line. First, let's go through adding a new file for your lab notebook. You can use your favorite text editor to make a new file, but all examples here will use `emacs` (because of course it is way better than vim, and definitely way cooler than nano. But if you already have a favorite editor then go for it).

```
emacs CSE185_Week2_LabNotebook.md
```
will create the file. Add some text (e.g. a title and the date). Use `ctrl-x ctrl-s` to save the file (`ctrl` refers to the control key). Now use `ctrl-x ctrl-c` to exit.

Now, let's add that file to the repository. This tells `git` to keep track of that file. You can also keep untracked files in this directory. Git won't pay attention to them unless they're added.
```
git add CSE185_Week2_LabNotebook.md
```

Type `git status` from anywhere in the repository folder to see what changes have happened since the last commit. You should see that you have added a new file (this will likely be in green on your terminal). To tell Git to update our local repository with this change, we need to commit it, along with a short message describing the commit:
```
git commit -a -m"Adding lab notebook for Week 2"
```
The `-m` option is to give a descriptive message about your commit. The `-a` option should be used pretty much any time you commit, since it tells `git` to automatically include any modifications to the repository you've made in the commit. Now if you type 'git status' you will see a message that there is nothing to commit.

Finally, to make these changes actually visible in the repository on Github, we have to "push" them:
```
git push
```
This will ask for your Github username and password. After you push, go back to the web browser and see if your changes appear there. You should be able to see (in blue at the top of the list of files) when your last commit was. You can click on that to see exactly what changed.

If you ever make changes on the web browser and need to update the clone of your repository on `ieng6.ucsd.edu`, you can "pull" those changes using:

```
git pull
```

You are encouraged to edit your files on the command line, although we won't know the difference if you edit on the browser vs. the command line. Next week we'll write our own script that will need to be included in the submission. It will be pretty annoying to edit that from the web browswer so it's worth getting comfortable with the command line git workflow this week. 

**Reminder: unless you do "git push" we will not see your changes. So make sure you commit and push all changes before the deadline on Tuesday! Whatever you see in the web browser for your repository is what we will see **

<blockquote>
**UNIX TIP**: It is annoying to keep opening and closing a text file to make edits. One solution is to have two open terminal screens, one where your lab notebook is open for editing and another where you are running commands. Alternatively, you can do all this is one terminal screen by taking advantage of the UNIX concepts of "foreground" and "background". If you are editing a file in emacs, type "ctrl-Z" to put emacs in the background. This will keep it running, but return you to the terminal where you can enter new commands. Once you've run your next command and would like to copy it back to emacs, type "fg" which continues the most recently stopped job by bringing it back to the foreground. This will take you back to emacs where you can keep editing where you left off. Using "ctrl-Z" and "fg" allow you to easily copy and paste commands back and forth from your notebook to the terminal and vice-versa.
</blockquote>

## 2. Inspect the data from your roommate

This sequencing data is in the `public/week2` directory. Your roommate’s data is labeled
`roommate.fastq`. Record how many reads there are in this file, then look at the first 20 lines with the
head command and answer the IClicker question.

Reminder - to get from your home directory to the parent cs185s class directory:
```shell
cd ..
```

Since there seem to be reads of multiple lengths, you suspect the data may have been processed
somehow. You can use a bash one-liner with `awk`, `sort`, and `uniq` to find out. Below is a reminder for
the general syntax of awk.

<blockquote>
**UNIX TIP**: A pipe (`|`) after a command redirects the output of the previous command as input to the next command. 
</blockquote>

```shell
awk '/search_pattern/ {actiontotakeonmatches; otheractions;}' file_to_awk
```

For the actual command (below), use pipes (`|`) to send the output of one command as input to the next. 

```shell
cat roommate.fastq | awk 'NR%4==0 {print length}' | sort -n | uniq -c
```

Let's walk through the commands we piped together:
* `cat`: This command, which stands for "concatenate", just outputs the contents of a file. Here we are opening the SRR file with `cat` and piping it into `awk`, so we don’t have to specify the file to `awk` after the brackets. 
* `awk`: `NR%4`, means NR ("number of records", or line number) modulo 4. The command returns every 4th line (because 4 goes into 4, 8, 12, 16, etc. perfectly with no remainder). This gives us a quick way to scan through each read (remember each read has 4 lines). The 4th line is the quality string, with a value for each base, so its length is the same as the length of the read.
* `sort` sorts the output lines. The `-n` flag tells it it's dealing with numbers so sort numerically. The next command requires the input to be sorted first.
* `uniq` returns only unique lengths (many reads are probably the same length), and the `-c` flag gives a count of how many there are of each.

**Record** the maximum read length, then answer the IClicker again. It looks like the reads come in various lengths, which suggests that someone has already processed the data. You check
with your friends in the med school, and it turns out that yes, they trimmed the low quality bases from
the ends of the reads for you.

## 3. Align your roommate's data to the reference sequence

The reference sequence for the influenza hemagglutinin gene is not in the public folder, so you will have to
download it from NCBI. This can be done with a command from the EntrezDirect utility, which is
installed for this course. 

`cd` into the `week2` directory you created above, and use the EntrezDirect command `efetch` to download the reference sequence, which has the NCBI id number: KF848938.1. An example of the `efetch` usage is below, you can also type `efetch -help` for more information. You must specify which NCBI database to use, the id or accession number of the sequence, and the format you would like. Redirect the output into a file with the “.fasta” extension.

<blockquote>
**UNIX TIP**: Using the symbol `>` after a command redirects the output (called "standard output") to a file rather directly to the screen.
</blockquote>

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
bwa mem KF848938.1.fasta /pathto/roommate.fastq | \
    samtools view -S -b | \
    samtools sort > roommate.bam
```

Run samtools view with the `-f 4` command to extract all unmapped reads, then count them. Use the number of unmapped reads and the total reads you calculated from the unprocessed `roommate.fastq` file to **calculate the number of reads that mapped.** 

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

Run VarScan on the mpileup. You can use the `VarScan` jar file in the `public/bin` directory. First, look for positions where most of the viruses infecting your roommate differ from the reference. Maybe there is a common mutation that wouldn’t have shown up in the HI test. Use a high minimum variant frequency cut-off to find only those mutants present in most (95% or more = 0.95) of the viral DNA molecules. 

```shell
java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp \
    roommate.mpileup \
    --min-var-freq 0.95 \
    --variants --output-vcf 1 > roommate.vcf
```

How many variants are reported back? To clean up the data and allow you to easily copy and paste,
**edit the awk script below so it only outputs the fields you are interested in** (position, reference base,
alternate base, in that order). Run the command, then **record the reference base, position, and mutant
base in your lab notebook.**

<blockquote>
**UNIX TIP**: The grep tool allows you to easily filter lines using [regular expressions](http://www.rexegg.com/regex-quickstart.html). The grep command below means to filter all lines starting with a "#" symbol (header lines).
</blockquote>

```shell
cat roommate.vcf | grep -v "^#" | {print $1, $2}' # Note you'll need to edit this
```

What do these mutations do? Could they be what allowed your roommate’s virus to escape the antibodies in your body from the flu vaccine? Since we are only looking at a single gene (and since the VEP doesn’t have a function for viruses) we will do this manually. In a browser, open the free, online sequence editor called WebDSV at
http://www.molbiotools.com/WebDSV/index.html

Go to NCBI (in a browser) and search for the reference sequence KF848938.1. Click on the FASTA link in the first result, then copy and paste the DNA sequence (not the header) where it says ‘paste sequence here’ on WebDSV. Alternatively you can copy the sequence from the reference fasta file you already downloaded to your directory. Click the yellow "process" button, then click "select all" and "translate". Close the pop up window, and when you get back to WebDSV, you should see the single-letter amino acid code above your sequence in the bottom window.

In the right side of the app, click the button with "AA" on it to open a codon table. For each position in your list, find the codon it is a part of. (You can hover over the bases to show the position, and you can hover over the amino acids to show the codon). Use the codon table to translate each mutation, and **record the original codon, mutated codon, original amino acid, its position in the protein, and the mutated amino acid. Then, record
whether the change is synonymous or missense.**

Example:
A72G ACA>ACG Thr24Thr synonymous

Repeat with the rest of the variants, record the results, and answer the IClicker question.

## 6. Look for rare variants with VarScan

Now try looking for rare variants. Set the min var freq to 0.001 (0.1%) and run the scan again, on the
same mpileup file. 

```shell
java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp \
    roommate.mpileup \
    --min-var-freq 0.001 \
    --variants --output-vcf 1 > roommate_rare.vcf
```

How many variants are reported back now, and how abundant are they? Modify the `awk` script again
to extract the position, reference base, and frequency. The frequency is part of field `$10`, so run the
command below.

```shell
cat roommate_rare.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}'
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
