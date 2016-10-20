import gzip

#opens the files
file=open('ex.fastq','r')
output=open('goodseq.fastq.gz','w')
not_use=open('badseq.fastq.gz','w')

#Finds the number of lines, as well as create a matrix of all the lines.
num_lines = sum(1 for line in open('ex.fastq'))
lines=file.readlines()
prob_list=[]

#Creates a list called prob_list of all the probabilities of a correct 18bp barcode.
for i in range(1,num_lines/4+1):
    score_string=lines[4*i-1][0:18]
    score_string = score_string.strip()
    prob=1
    for char in score_string:
        prob=prob*(1-(10**(-(ord(char)-33)/10)))
    prob_list.append(str(prob))

#Creates the lines to be added to the new file with only index, barcode, and score.
def update_file(ln, score,lines,output):
    for i in range(0,len(lines)):
        if i == ln:
            index=lines[i-1][-7:].strip('\n')
            new_line = lines[i][0:18].strip('\n') + score +'\n'
    return index, new_line

#Iterates through the lines and write one file called good seq with all seqs with scores > 0.9
#and all the other seqs to another file called badseq.
for i in range(1,(num_lines/4)+1):
    index,new_line=update_file(4*i-3,prob_list[i-1],lines,output)
    if float(prob_list[i-1]) > 0.9:
        output.write('{}\n{}'.format(index,new_line))
    else:
        not_use.write('{}\n{}'.format(index,new_line))
