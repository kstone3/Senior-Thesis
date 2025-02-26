#!/usr/bin/env python
"""
Simulate the evolution of South Cascade Glacier (RGI60-11.00897)
from 1984 to 2024 using OGGM.
"""

import os
import matplotlib.pyplot as plt

import oggm
from oggm import cfg, workflow, tasks

# ---------------------------
# 1. OGGM Setup and Configuration
# ---------------------------

# Initialize OGGM (set logging level as desired)
oggm.cfg.initialize(logging_level='INFO')

# Define a working directory (adjust path as needed)
workdir = '../OGGM'
if not os.path.exists(workdir):
    os.makedirs(workdir)
oggm.cfg.PATHS['working_dir'] = workdir

# Optionally, disable multiprocessing for reproducibility
cfg.PARAMS['use_multiprocessing'] = False

# ---------------------------
# 2. Glacier Initialization
# ---------------------------
# South Cascade Glacier is identified here by its RGI ID.
# (This RGI id is used in several OGGM examples.)
rgi_id = 'RGI60-02.18778'

# Initialize the glacier region (glacier directory) for the selected glacier
gdir = workflow.init_glacier_region(rgi_id=rgi_id, path=workdir)

# ---------------------------
# 3. Pre-Processing Tasks
# ---------------------------
# Run the standard processing tasks to compute centerlines, catchments, etc.
workflow.execute_entity_tasks([gdir], [
    tasks.define_glacier_region,
    tasks.glacier_masks,
    tasks.compute_centerlines,
    tasks.catchment_area,
    tasks.catchment_intersections,
    tasks.compute_downstream_line,
    tasks.compute_inversion_altitude,
    tasks.compute_ref_t_star,
    tasks.compute_ref_t_bias  # if bias correction is desired
])

# ---------------------------
# 4. Glacier Evolution Model
# ---------------------------
# Run the glacier evolution model by calling the `run_glacier_flowline` task
# The simulation will be run from 1984 to 2024, using available climate data.
tasks.run_glacier_flowline(gdir, start_year=1984, end_year=2024)

# ---------------------------
# 5. Read and Plot Model Output
# ---------------------------
# The model run saves its output (typically a dict with keys like 'time', 'volume', etc.)
model_output = gdir.read_pickle('model_run_output')

# Extract time series and volume evolution
times = model_output['time']
volumes = model_output['volume']

plt.figure(figsize=(10, 6))
plt.plot(times, volumes, label='Glacier Volume')
plt.xlabel('Year')
plt.ylabel('Volume (m³)')
plt.title('Evolution of South Cascade Glacier Volume (1984-2024)')
plt.legend()
plt.grid(True)
plt.show()

print("Final glacier volume in 2024: {:.2f} m³".format(volumes[-1]))