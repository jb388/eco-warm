This directory contains template tables for the Synthesizing Warming Experiments to Depth Data Integration Effort (SWEDDIE).

These templates correspond to the "core" component of the database, i.e., data that is predominantly static. The data tables should be self explanatory; the accompanying 'xxxMeta' tables define the variables for each table. Briefly, the 'experiment' table focuses on infrastructure; 'site' focuses on site characteristics, e.g., long term climate averages, and soil, vegetation, and geological characterization; 'plot' contains experimental design information.
Update 26-May-2024: change 'xxxMeta' naming convention to '_dd' to match CSV specifications from ESS-DIVE (this is cosmetic, but may help avoid confusion about metadata versus data dictionary/annotation files)

We have also provided template tables ('datHeatTS_template', 'datHeatTS_ann_template') for reporting timeseries of active heating. The tables are designed to be filled out per plot, i.e., each plot will have a separate 'datHeatTS' table, with the accompanying metadata/annotation 'datHeatTS_ann' table describing the variables. If you use the template files, please replace the "_template" component of the table name with a plot ID, e.g., 'datHeatTS_plot1', 'datHeatTS_ann_plot1'. Note that if you have this data in another format, or if data for all the plots are in the same file, it is fine to submit as is, as long as the submitted file contains the dates for which the heating infrastructure was active for each plot, and the plot IDs correspond to those entered in the "plot" table.

Notes:
1) measurement data (active heating timeseries, temperature, moisture, CO2 fluxes, etc.) are linked to core tables via the 'exp_name', 'sit_name', and 'plt_name' variables
2) for multiple entries within a single cell, please separate with ";" (semicolon) and no spaces
3) Missing data values preferably "NA". Please specify missing data values, not a number values (e.g., "NAN"), etc., in the "dat" metadata files/annotation files

More information about the database is available here.


Please contact Jeff Beem-Miller (jbeemmil@umich.edu) with any questions.

Thanks!
