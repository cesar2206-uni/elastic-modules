# Libraries
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from io import StringIO

# Matplotlib configuration
mpl.rcParams['font.family'] = 'Yu Gothic'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2


# Sidebar Layout
st.sidebar.markdown("# Upload Data")

uploaded_file = st.sidebar.file_uploader("Choose a .csv file:")

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file, names = ["strain", "stress"])
    #string_value = StringIO(uploaded_file.getvalue().decode("utf-8")).read()
correction_0 = st.sidebar.checkbox('Correction 1: XY in (0, 0)')
correction_mpa = st.sidebar.checkbox('Correction 2: Y axis in mPa')
st.sidebar.markdown("# Secant elastic modulus")

correction_secant = st.sidebar.checkbox('Correction of Secant elastic modulus')
d_secant = st.sidebar.slider('Initial displacement', 0.0, 2.0, 1.0)

st.sidebar.markdown("# Linear elastic modulus")
fit_value = st.sidebar.number_input("Minimum fit lenght", min_value = 3, value = 7, step = 1)



# Processing
   
strain = df["strain"].to_numpy()
stress = df["stress"].to_numpy()

if correction_0:
    strain = np.concatenate((np.array([0]), df["strain"].to_numpy()))
    stress = np.concatenate((np.array([0]), df["stress"].to_numpy()))

def plot_modulus(fit_value):
    
    # Algorithm of chi square
    def chi_algorithm(x, y, fit_value):
        num_points = len(x)
        min_fit_length = fit_value
        chi = 0
        chi_min = 10000000
        i_best = 0
        j_best = 0

        for i in range(len(x) - min_fit_length):
            for j in range(i+min_fit_length, len(x)):

                coefs = np.polyfit(x[i:j],y[i:j],1)
                y_linear = x * coefs[0] + coefs[1]
                chi = 0
                for k in range(i,j):
                    chi += ( y_linear[k] - y[k])**2

                if chi < chi_min:
                    i_best = i
                    j_best = j
                    chi_min = chi
                    #print(chi_min)

        coefs = np.polyfit(x[i_best:j_best],y[i_best:j_best],1)
        y_linear = x[i_best:j_best] * coefs[0] + coefs[1]
        return x[i_best:j_best], y_linear

    # Calculation of tangent modulus
    EX_tan, EY_tan = chi_algorithm(strain, stress, fit_value)
    if correction_mpa:
        E_tan = (EY_tan[1] - EY_tan[0])/(EX_tan[1] - EX_tan[0]) * 100
    else:
        E_tan = (EY_tan[1] - EY_tan[0])/(EX_tan[1] - EX_tan[0]) * 100/1000
    # Obtain the secant stiffness modulus
    max_stress = max(stress)
    med_stress = max_stress / 2

    for i, stress_value in enumerate(stress):
        if med_stress >= stress[i] and med_stress <= stress[i + 1]:
            med_strain = strain[i] + (med_stress - stress[i]) * (strain[i + 1] - strain[i]) / (stress[i + 1] - stress[i])
            break
    
    if correction_secant:     
        for i, strain_value in enumerate(strain):
            if d_secant >= strain[i] and d_secant <= strain[i + 1]:
                d_stress = stress[i] + (d_secant - strain[i]) * (stress[i + 1] - stress[i]) / (strain[i + 1] - strain[i])
                break

    # Calculation of secant modulus
    
    if correction_secant:
        EX_sec = np.array([d_secant, med_strain])
        EY_sec = np.array([d_stress, med_stress])
        if correction_mpa:
            E_sec = (EY_sec[1] - EY_sec[0])/(EX_sec[1] - EX_sec[0]) * 100
        else:
            E_sec = (EY_sec[1] - EY_sec[0])/(EX_sec[1] - EX_sec[0]) * 100/1000    
    else:
        EX_sec = np.array([0, med_strain])
        EY_sec = np.array([0, med_stress])
        if correction_mpa:
            E_sec = (EY_sec[1] - EY_sec[0])/(EX_sec[1] - EX_sec[0]) * 100
        else:
            E_sec = (EY_sec[1] - EY_sec[0])/(EX_sec[1] - EX_sec[0]) * 100/1000

    # Plotting zone
    
    em_fig = plt.figure(figsize=(10, 8))
    em_ax = em_fig.add_subplot(111)

    em_ax.plot(strain, stress, "blue")
    em_ax.plot(EX_tan, EY_tan, color = "red")
    
    em_ax.text( 
        (EX_tan[0] + EX_tan[-1])/2 -max(strain)/10 , 
        (EY_tan[0] + EY_tan[-1])/2, 
        r'$E_{lin}$', fontsize=15, 
        color = "red"
        )

    em_ax.plot(EX_sec, EY_sec, "orange")
    if correction_secant:
        em_ax.plot([EX_sec[0], EX_sec[0]], [0, EY_sec[0]], linestyle = "dashed", color = "gray")
    #em_ax.plot([0, EX_sec[0]], [EY_sec[0], EY_sec[0]], linestyle = "dashed", color = "gray")
    em_ax.text(
        (EX_sec[0] + EX_sec[-1])/2 - max(strain)/10,
        (EY_sec[0] + EY_sec[-1])/2, 
        r'$E_{sec}$', fontsize=15,
        color = "orange"
        )


    em_ax.text(
        0.6 * max(strain), 
        0.1 * max(stress), 
        "Resultados: \n E. Lineal = %.1f MPa \n E. Secante = %.1f MPa" % (E_tan,E_sec),
        fontsize= 13, 
        bbox=dict(
            facecolor='white', 
            alpha=0.7,)
            )

    em_ax.xaxis.set_tick_params(
        which='major', 
        size=10, 
        width=2, 
        direction='in', 
        right='on'
        )

    em_ax.yaxis.set_tick_params(
        which='major', 
        size=10, 
        width=2, 
        direction='in', 
        right='on'
        )

    em_fig.suptitle('Ensayo de compresión no confinada')
    em_ax.set_title(uploaded_file.name[:-4])
    em_ax.set_xlabel("Deformación Axial, %", labelpad=10)
    
    if correction_mpa:
        em_ax.set_ylabel('Esfuerzo a compresión, mPa', labelpad=10)
    else:
        em_ax.set_ylabel('Esfuerzo a compresión, kPa', labelpad=10)
    
    return em_fig


# Principal Layout
st.markdown("# Evaluation of Elastic Modulus")   

#if correction_0:
em_result = plot_modulus(fit_value)
st.pyplot(em_result)

fn = uploaded_file.name[:-3]+ "png"
em_result.savefig(fn) 
with open(fn, "rb") as img:
    btn = st.download_button(
        label="Download image",
        data=img,
        file_name=fn,
        mime="image/png"
    )
        
st.markdown("Author: César Manuel Sánchez Oré (cesar.sanchez@wsp.com)")