#!/usr/bin/env python
# coding: utf-8

# In[1]:


# individual information add & rearrange
import sys

ii = sys.argv[1]

for i in range(int(ii),int(ii)+1):
    print 'chr' + str(i).replace("23","X").replace("24","Y")
    input_file = file("revision_CpG_CtoT_PCAWG_CGI_reclassification_chr%s.txt" %str(i).replace("23","X").replace("24","Y"))
    info_file = file("revision_CpG_CtoT_PCAWG_chr%s_individual.txt" %str(i).replace("23","X").replace("24","Y"))
    output_file = file("revision_CpG_CtoT_PCAWG_CGI_reclassification_chr%s_individual.txt" %str(i).replace("23","X").replace("24","Y"),"w")
    
    input_line = input_file.readline().strip()
    info_line = info_file.readline().strip()
    
    while input_line:
        info_split = info_line.split('\t')
        output_file.write(input_line + '\t' + '\t'.join(info_split[25:]) + '\n')

        input_line = input_file.readline().strip()
        info_line = info_file.readline().strip()    

    
input_file.close()

output_file.close()
print 'THE END'



# In[ ]:




