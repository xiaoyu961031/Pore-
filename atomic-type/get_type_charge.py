# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 16:15:51 2022

@author: Xiaoyu.Wu
"""

dir = r'qmof_cifs_with_charges/qmof_cifs_with_charges'

dict_elements = {
                  'B' :[1,1, [0],         ], # Dict {'Element': [group_numbers, sa_index, threshold]}
                  'Br':[2,2, [-0.3],      ],
                  'C' :[4,3, [-0.4,0,0.4],],
                  'Cl':[3,4, [-0.3,0],    ],
                  'F' :[2,5, [-0.3],      ],
                  'H' :[2,6, [0.2],       ],
                  'I' :[3,7, [-0.3,0],    ],
                  'N' :[4,8, [-0.5,0,0.3],],
                  'O' :[3,9, [-0.8,-0.4], ],
                  'P' :[3,10,[0.6,1.2],   ],
                  'S' :[4,11,[0,0.5,1.0], ],
                  'Si':[2,12,[1,2],       ],
                }

MOFlist_clean=[]
MOFlist = open('MOFlist.txt','r')
MOFlist_content = MOFlist.readlines()
for i in MOFlist_content:
    MOFlist_clean.append(i.strip('\r\n').strip('\n') )
    
element_SA_file = open('element_SA','r')
element_SA_file_content = element_SA_file.readlines()
element_type_file = open('element_atomic_SA','w')


def get_type_wise(list, sa, group_num, threshold_list):
    
    if group_num == 1:
        
        SA_list = [0]
        return SA_list  
    
    if group_num == 2:
    
        if sa != 0:
            threshold = threshold_list[0]
            
            class1 = []
            class2 = []
            
            for i in list:
                if float(i)<threshold:
                    class1.append(i)
                if float(i)>=threshold:
                    class2.append(i)
            
            SA_1 = float(sa)*(int(len(class1))/len(list))
            SA_2 = float(sa)*(int(len(class2))/len(list))
            
            SA_list = [str(round(SA_1, 3)), str(round(SA_2, 3))]
            return SA_list        
        
        else:
            SA_list = ['0','0']
            return SA_list

    if group_num == 3:    
    
        if sa != 0:    
            threhold1 =  threshold_list[0]
            threhold2 =  threshold_list[1]
        
            class1 = []
            class2 = []
            class3 = []
            
            for i in list:
                if float(i)<threhold1:
                    class1.append(i)
                if float(i)>= threhold1 and float(i)<threhold2:
                    class2.append(i)
                if float(i)> threhold2:
                    class3.append(i)                      
            
            SA_1 = float(sa)*(int(len(class1))/len(list))
            SA_2 = float(sa)*(int(len(class2))/len(list))
            SA_3 = float(sa)*(int(len(class3))/len(list))
            
            SA_list = [str(round(SA_1, 3)), str(round(SA_2, 3)), str(round(SA_3, 3))]
            return SA_list        
        
        else:
            SA_list = ['0','0','0']
            return SA_list

    if group_num == 4:             
    
        if sa != 0:    
            threhold1 =  threshold_list[0]
            threhold2 =  threshold_list[1]
            threhold3 =  threshold_list[2]
        
            class1 = []
            class2 = []
            class3 = []
            class4 = []
            
            for i in list:
                if float(i)<threhold1:
                    class1.append(i)
                if float(i)>= threhold1 and float(i)<threhold2:
                    class2.append(i)
                if float(i)>= threhold2 and float(i)<threhold3:
                    class3.append(i)      
                if float(i)>= threhold3:
                    class4.append(i)                 
            
            SA_1 = float(sa)*(int(len(class1))/len(list))
            SA_2 = float(sa)*(int(len(class2))/len(list))
            SA_3 = float(sa)*(int(len(class3))/len(list))
            SA_4 = float(sa)*(int(len(class4))/len(list))
            
            SA_list = [str(round(SA_1, 3)), str(round(SA_2, 3)), str(round(SA_3, 3)), str(round(SA_4, 3))]
            return SA_list
        
        else:
            SA_list = ['0','0','0','0']
            return SA_list
        
    else:  
        print ('Please mannually set up group '+str(group_num)+' classification!'+'\n')


for mof in MOFlist_clean:
        
    cif_path = dir+'/'+mof+'.cif'
    cif =  open(cif_path,'r')
    cif_content = cif.readlines() 
    
    print (str(mof)+' ')
    element_type_file.writelines(str(mof)+' ')
    
    for cif_line in cif_content:
        if '_atom_site_charge' in cif_line:
            tmpindex = cif_content.index(cif_line)
            atomic = cif_content[tmpindex+1:]
    
    for key,value in dict_elements.items():

        element_type = str(key)
        group_num = value[0]
        sa_index = value[1]
        threshold_list = value[2]
        element_charge_list = []
            
        for i in atomic:
            if element_type+' ' in i:                
                element_charge = i.split()[7]
                element_charge_list.append(element_charge)  
                    
        for i in element_SA_file_content:
            if str(mof) == str(i.split()[0]):
                element_sa = float(i.split()[sa_index])
            
        sa_type_list = get_type_wise(element_charge_list, element_sa, group_num, threshold_list)
                
        for i in sa_type_list:
            print (str(i)+' ')
            element_type_file.writelines(str(i)+' ')

        
    print ('\n')
    element_type_file.writelines('\n')         
            
        
           
        
            


    
    
    


