#This little script it is useful for search and recover the FimH sequence in a database of FASTA CDS sequences, little modificatiuon will be allow to search any gen.
#This script must be located in the directory with the files of the sequence in a .txt format.

#Create a folder called "output"
mkdir ./output

#Create a tect file called "results"

touch ./output/FimH-bp.txt

#For a variable "i" in the name of the text file inside of the directory

for i in *.txt;

#Print the name of the file, grep find "FimH" and print 12 lines. 

do echo ">"$i >> ./output/resultados.txt
grep FimH -A13 -h $i | tail -n +2  >> ./output/resultados.txt

done

echo "done"
