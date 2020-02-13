#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# thousand + PCAWG analysis for individual variation (Except Y chromosome)

list_class = ['noCGI', 'NR_TSS', 'NR_intragenic', 'intergenic', 'NM_TSS', 'NM_intragenic', 'miscellaneous']
total_class = [0,0,0,0,0,0,0]

individual_class = []
for k in range(30,4357):
    individual_class.append([0,0,0,0,0,0,0])
    
for i in range(1,24):
    print 'chr' + str(i).replace("23","X").replace("24","Y")
    input_file = file("revision_CpG_CtoT_PCAWG_CGI_reclassification_chr%s_individual.txt" % str(i).replace("23","X").replace("24","Y"))
    input_line = input_file.readline().strip()
    input_line = input_file.readline().strip()
    while input_line:
        input_split = input_line.split('\t')
        input_gt_pcawg = input_split[24].split(';')[0]
        input_gt_thousand = input_split[24].split(';')[1]
        
        if input_split[25] == 'x':
            input_class = 'noCGI'
        else:
            input_class = input_split[29]
            if len(input_class.split(',')) > 1:
                input_class = 'miscellaneous' 
                         
        c_index = list_class.index(input_class)
        
        total_class[c_index]+=1
        for j in range(30,1853): #PCAWG
            if input_split[j] == '.':
                'blank'
            else:
                for kk in input_split[j].split('|'):
                    if kk == input_gt_pcawg:
                        individual_class[j-30][c_index] += 1
                        break
                    
            
        for j in range(1853,4357): #Thousand
            if input_split[j] == '.':
                'blank'
            else:            
                for kk in input_split[j].split('|'):
                    if kk == input_gt_thousand:
                        individual_class[j-30][c_index] += 1
                        break

                
        input_line = input_file.readline().strip()
    
print list_class
print total_class
print individual_class[0:3]


output_file = file("figure1e_individual_union.txt","w")
output_file.write('\t'.join(list_class) + '\n')
output_file.write('\t'.join(map(str,total_class)) + '\n')
for k in range(30,4357):
    output_file.write('\t'.join(map(str,individual_class[k-30])) + '\n')

output_file.close()

print 'THE END'    

