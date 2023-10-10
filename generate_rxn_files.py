from lxml import etree
#import argparse
import sys


file_tg =  "Rxn_module_TG.xml"
file_caoutbuf = "Rxn_module_CaOutBuf.xml"
file_cabuf = "Rxn_module_CaBuf.xml"
file_fluo = "Rxn_module_Fluo4.xml"
file_fura ="Rxn_module_Fura2.xml"

file_list_full_ER_buf_SERCA = [                 
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_mGLuR.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    "Rxn_module_TG.xml"
]
file_list_ER_no_mGLur_buf_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
file_list_ER_no_mGLur_no_RyR_buf_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
file_list_ER_no_mGLur_no_RyR_no_SOCE_buf_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
]
#"For checking if adding SOCE help Ca wave propagation
file_list_ER_no_IP3R_buf_SERCA = [
    "Rxn_module_Ca.xml",
    #    "Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    ]

#for reproducing Ca wave paper
file_list_ER_no_IP3R_no_SOCE_buf_SERCA = [
    "Rxn_module_Ca.xml",
    #    "Rxn_module_CaBuf.xml",
    "Rxn_module_Fura2.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_buffering_SERCA_ER.xml",
    ]

file_list_full_ER_simple_SERCA = [                 
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_mGLuR.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    "Rxn_module_TG.xml"
]
file_list_ER_no_mGLur_simple_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
file_list_ER_no_mGLur_no_RyR_simple_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
file_list_ER_no_mGLur_no_RyR_no_SOCE_simple_SERCA = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_simple_SERCA_ER.xml",
]
#"For checking if adding SOCE help Ca wave propagation
file_list_ER_no_IP3R_simple_SERCA = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    ]

#for reproducing Ca wave paper
file_list_ER_no_IP3R_no_SOCE_simple_SERCA_Fura2 = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_Fura2.xml",
    ]

file_list_ER_no_IP3R_simple_SERCA_Fura2 = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    "Rxn_module_Fura2.xml",

    ]

#for reproducing Ca wave paper
file_list_ER_RyR2CaM_no_SOCE_simple_SERCA_Fura2 = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_modified.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_Fura2.xml",
    ]


#"For checking if adding SOCE help Ca wave propagation
file_list_ER_RyR2CaM_simple_SERCA = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_CaM.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    ]

#for reproducing Ca wave paper
file_list_ER_RyR2CaM_no_SOCE_simple_SERCA = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_CaM.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    ]

file_list_ER_RyR2CaM_simple_SERCA_Fura2 = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_CaM.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    "Rxn_module_Fura2.xml",

    ]

#for reproducing Ca wave paper
file_list_ER_RyR2CaM_no_SOCE_simple_SERCA_Fura2 = [
    "Rxn_module_Ca.xml",
    #"Rxn_module_CaBuf.xml",
    "Rxn_module_RyR_KeizerSmith_CaM.xml",
    "Rxn_module_simple_SERCA_ER.xml",
    "Rxn_module_Fura2.xml",
    ]


def read_in_files(flist):
    roots = []
    species = set()
    specie_diff = {}
    rxn_ids = set()
    my_rxn_file = etree.Element("ReactionScheme")

    for fname_xml in flist:
        print(fname_xml)
        tree = etree.parse(fname_xml)
        roots.append(tree.getroot())
        for son in roots[-1]:
            if son.tag == "Specie":
                specie_id = son.attrib["id"]
                specie_kdiff = son.attrib["kdiff"]
                species.add(specie_id)
                if specie_id in specie_diff:
                    if specie_diff[specie_id]!= specie_kdiff:
                        print(specie_id, specie_diff[specie_id], specie_kdiff)
                else:
                    specie_diff[specie_id] = specie_kdiff
                    my_rxn_file.append(son)
            elif son.tag == "Reaction":
                if son.attrib["id"] in rxn_ids:
                    print(son.attrib["id"], " already exists", fname_xml)
                else:
                    rxn_ids.add(son.attrib["id"])
                    my_rxn_file.append(son)
    return my_rxn_file

if __name__ == "__main__":
    
    # 1 no mGluR no RyR
   

    my_rxn_f = read_in_files(file_list_ER_RyR2CaM_simple_SERCA_Fura2)
    with  open("Rxn_RyR2CaM_SERCA_simple_Fura2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    
    my_rxn_f = read_in_files(file_list_ER_RyR2CaM_no_SOCE_simple_SERCA_Fura2)
    with open("Rxn_RyR2CaM_no_SOCE_SERCA_simple_Fura2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))

    my_rxn_f = read_in_files(file_list_ER_RyR2CaM_simple_SERCA_Fura2)
    with open("Rxn_RyR2CaM_SERCA_simple.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(file_list_ER_RyR2CaM_no_SOCE_simple_SERCA_Fura2)
    with open("Rxn_RyR2CaM_no_SOCE_SERCA_simple.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f,
                               pretty_print=True).decode("utf-8"))
                    

   
    my_rxn_f = read_in_files(file_list_ER_no_IP3R_simple_SERCA_Fura2)
    with open("Rxn_no_IP3R_SERCA_simple_RyR2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))

    my_rxn_f = read_in_files(file_list_ER_no_IP3R_simple_SERCA_Fura2)
    with open("Rxn_no_IP3R_SERCA_simple_RyR2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))


    my_rxn_f = read_in_files(file_list_ER_no_IP3R_no_SOCE_simple_SERCA_Fura2)
    with open("Rxn_no_IP3R_no_SOCE_SERCA_simple_RyR2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))

    my_rxn_f = read_in_files(file_list_ER_no_IP3R_no_SOCE_simple_SERCA_Fura2)
    with open("Rxn_no_IP3R_no_SOCE_SERCA_simple_RyR2_Fura2.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
