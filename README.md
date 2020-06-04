# FFD_PFM_test_bed
A test-bed for FFD (Fast Fluid Dynamics) and PFM (Potential Flow Model) methods for simulating velocity and pressures

A brief summary of the files in this project is as follows:

1. main.py
This contains the GUI of the test bed. It is connected to the solvers - CFD_solver.py and PFM_solver.py

2. FFD_solver.py
This contains the Fast Fluid Dynamics model to compute the velocity and pressures in a staggered grid structure.

3. PFM_solver.py
This contains the Potential Flow Model to compute the velocity and pressure values from a potential flow model.

## Running the Test bed and visualizing the results

### Python Dependencies
* `python 3.7`
* `dash 0.43`
* `numpy latest`

`pip install dash==0.43`
`pip install numpy`

### Running on Local host (Windows/Linux)
Run the following command in the current directory through command-line interface
`python main.py`

Or run the `main.py` script through python IDE.

Open following link on web browser:
http://127.0.0.1:8050/
