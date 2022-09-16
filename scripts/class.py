import streamlit as st

st.write("""
# Genotype Calling Application

## Chromatogram

## Peak range table


## Star allele calling
""")


allele_data = {
    'CYP2D6_14': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (25,35)
                        },
                    2: {    'type': 'mutant',
                            'base': 'A',
                            'range': (27,36)
                        } 
              },
    'CYP2D6_10B': { 1: {    'type': 'widtype',
                            'base': 'C',
                            'range': (28,38),
                        },
                    2: {    'type': 'mutant',
                            'base': 'T',
                            'range': (31,37),
                        } 
            },
    'CYP2D6_49': { 1: {    'type': 'widtype',
                            'base': 'T',
                            'range': (37,44),
                        },
                    2: {    'type': 'mutant',
                            'base': 'A',
                            'range': (38,45),
                        } 
            },
    'CYP2D6_21': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (41,45)
                        },
                    2: {    'type': 'mutant',
                            'base': 'C',
                            'range': (42,49),
                        } 
            },
    'CYP2D6_41': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (45,50),
                        },
                    2: {    'type': 'mutant',
                            'base': 'A',
                            'range': (47,50),
                        } 
              },
    'CYP2D6_52': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (50,55),
                        },
                    2: {    'type': 'mutant',
                            'base': 'A',
                            'range': (50,55),
                        } 
            },
    'CYP2D6_18': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (53,57),
                        },
                    2: {    'type': 'mutant',
                            'base': 'T',
                            'range': (55,59),
                        } 
            },
    'CYP2D6_2': { 1: {    'type': 'widtype',
                            'base': 'C',
                            'range': (57,62),
                        },
                    2: {    'type': 'mutant',
                            'base': 'T',
                            'range': (59,63),
                        } 
            },
    'CYP2D6_60': { 1: {    'type': 'widtype',
                            'base': 'G',
                            'range': (63,68),
                        },
                    2: {    'type': 'mutant',
                            'base': 'A',
                            'range': (65,68),
                        } 
            },
    'CYP2D6_5': { 1: {    'type': 'widtype',
                            'base': 'A',
                            'range': (68,73),
                        },
                    2: {    'type': 'mutant',
                            'base': 'G',
                            'range': (68,72),
                        } 
            },
}

def generate_allele_data(marker):
    data = allele_data.get(marker)
    
    for i, row in enumerate(data, start=1):
        row_data = data.get(row)
        type_base = f"{row_data['base']} ({row_data['type']})"
        st.write(f'Allele: {type_base}', color = 'blue')
        st.number_input("Min Bin", value = row_data['range'][0], key = f'{marker}_{i}_min')
        st.number_input("Max Bin", value = row_data['range'][1], key = f'{marker}_{i}_max')


with st.sidebar:
    upload = st.file_uploader('Upload FSA file', type='FSA', accept_multiple_files=True)
    st.button("Process")

    st.write('Adjust bin range')

    selected_marker = st.selectbox("Select Marker", allele_data.keys())
    generate_allele_data(selected_marker)




for file in upload:
    panel = file.name.split('-')[0]
    with st.expander(f'Panel {panel}'):
        st.write(file.name)

