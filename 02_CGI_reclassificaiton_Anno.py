#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys, os

ii = sys.argv[1]
info_file = file("CGI_classification_revised.txt")
info_line = info_file.readline().strip()
info_list = []
for i in range(1,25):
    info_list.append([])
    
while info_line:
    info_split = info_line.split('\t')
    info_chr = info_split[1]
    info_index = int(info_chr.replace("chr","").replace("X","23").replace("Y","24")) - 1
    info_start = int(info_split[3]) + 1
    info_end = int(info_split[4])
    info_CGI = info_chr + ':' + info_split[3]
    
    info_list[info_index].append([info_CGI,info_start,info_end,info_split[6],info_split[12]])
    info_line = info_file.readline().strip()

print 'info_list has been completed'
info_file.close()

for i in range(int(ii),int(ii)+1):
    print 'chr' + str(i).replace("23","X").replace("24","Y")
    input_file = file("revision_CpG_CtoT_PCAWG_chr%s.txt" %str(i).replace("23","X").replace("24","Y"))
    output_file = file("revision_CpG_CtoT_PCAWG_CGI_reclassification_chr%s.txt" %str(i).replace("23","X").replace("24","Y"),"w")
    
    input_line = input_file.readline().strip()
    output_file.write(input_line+'\tIn_CGI\tnear_CGI\tdistace_from_CGI\tCGI_size\tCGI_class\n')
    
    input_line = input_file.readline().strip()
    prev_CGI = ''
    while input_line:
        input_split = input_line.split('\t')
        input_chr = int(input_split[0].replace("chr","").replace("X","23").replace("Y","24"))
        input_pos = int(input_split[1])
        info = ['','','','','']
        
        input_distance = 999999999
        input_decision = ''
        for [info_CGI,info_start,info_end,info_size,info_class] in info_list[input_chr-1]:
            if input_pos >= info_start and input_pos <= info_end:
                info = ['o',info_CGI,'0',info_size,info_class]                
                break
            else:
                temp_distance1 = input_pos - info_start
                temp_distance2 = input_pos - info_end
                if min(abs(temp_distance1),abs(temp_distance2)) == abs(temp_distance1):
                    temp_distance = temp_distance1
                else:
                    temp_distance = temp_distance2
                
                if input_distance > abs(temp_distance):
                    input_distance = abs(temp_distance)
                    info = ['x',info_CGI,str(temp_distance),info_size,info_class]
        
        output_file.write(input_line+'\t' + '\t'.join(info) + '\n')
        input_line = input_file.readline().strip()
        
    
    
    input_file.close()
    output_file.close()

print 'THE END'

