import streamlit as st
import matplotlib.pyplot as plt
from Wave_Packet import wave_packet

#Streamlit App
st.title("Quantum Tunnelling with Gaussian Wavepacket")

# Sidebar
with st.sidebar:
    st.header("Simulation Parameters")
    sigma = st.slider("Initial width of the wavepacket", 0.1, 10.0, 5.0, 0.1)
    x0 = st.slider("Initial position of the wavepacket", -200, 200, -150, 1)
    k0 = st.slider("Initial wavevector of the wavepacket", 0.1, 10.0, 1.0, 0.1)
    x_begin = st.slider("Starting position of the x-axis", -500, 0, -200, 10)
    x_end = st.slider("Ending position of the x-axis", 0, 500, 200, 10)
    dt = st.slider("Time step", 0.01, 1.0, 0.01, 0.01)
    barrier_height = st.slider("Height of the potential barrier", 0.1, 10.0, 1.0, 0.1)
    barrier_width = st.slider("Width of the potential barrier", 0.1, 10.0, 4.0, 0.1)
    n_points = st.slider("Number of points on the x-axis", 100, 1000, 500, 100)
    num_steps = st.slider("Number of time steps", 100, 1000, 500, 100)

start_button = st.button("Start Simulation")
# Main content
wave_packet = wave_packet(sigma0=sigma, x0=x0, k0=k0, x_begin=x_begin, x_end=x_end, dt=dt, barrier_height=barrier_height, barrier_width=barrier_width, n_points=n_points)


# Display the animation
if start_button:
    fig,ax = plt.subplots()
    ax.plot(wave_packet.xgrid, wave_packet.potential*.1, color='r')
    ax.set_ylim(0, barrier_height*0.11)
    ax.set_xlim(x_begin, x_end)
    ax.set_ylabel('Probability density (a$_0$)')
    ax.set_xlabel('Position (a$_0$)')
    ax.grid(True)

    line, = ax.plot(wave_packet.xgrid, wave_packet.evolve(), label='Wavepacket')
    ax.legend()

    plot_placeholder = st.empty()

    #Run the simulation
    for _ in range(num_steps):
        data = wave_packet.evolve()
        line.set_ydata(data)
        
        #Redraw the plot
        fig.canvas.draw()
        plot_placeholder.pyplot(fig)




