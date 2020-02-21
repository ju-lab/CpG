#!/usr/bin/env python
# coding: utf-8

# In[4]:


import sys, os, gzip
chr_num = sys.argv[1]
#chr_num = '1'
output_file = file("union_signature_chr%s.txt" % chr_num.replace("23","X"),"w")

info_file = file("signature_1000_chr%s.txt" %chr_num.replace("23","X"))   
p_file = file("signature_pcawg_chr%s.txt" %chr_num.replace("23","X"))

info_line = info_file.readline().strip()
p_line = p_file.readline().strip()
info_line = info_file.readline().strip()
p_line = p_file.readline().strip()

output_file.write('#CHROM\tPOS\tID\tREF\tALT\tAC\n')

while info_line:
    info_split = info_line.split('\t')
    info_pos = int(info_split[1])
    
    if p_line != '':
        p_split = p_line.split('\t')
        p_pos = int(p_split[1])

        if info_pos == p_pos:
            output_file.write(info_line+'\n')
            output_file.write(p_line+'\n')
            info_line = info_file.readline().strip()
            p_line = p_file.readline().strip()
        elif info_pos < p_pos:
            output_file.write(info_line+'\n')
            info_line = info_file.readline().strip()
        elif p_pos < info_pos:
            output_file.write(p_line+'\n')
            p_line = p_file.readline().strip()
    else:
        output_file.write(info_line+'\n')
        info_line = info_file.readline().strip()        

while p_line:
    output_file.write(p_line+'\n')
    p_line = p_file.readline().strip()    

    
output_file.close()
print 'THE END'


# In[ ]:




