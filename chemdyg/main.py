import jinja2
import os
import pprint

from zppy.utils import (
    checkStatus,
    getComponent,
    getTasks,
    getYears,
    makeExecutable,
    submitScript,
)

# -----------------------------------------------------------------------------
def chemdyg(path, config, scriptDir, existing_bundles, job_ids_file):

    # Initialize jinja2 template engine
    path_extra = os.path.join(path,"templates")
    templateLoader = jinja2.FileSystemLoader(
        searchpath=(config["default"]["templateDir"], path_extra)
    )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # --- List of chemdyg tasks ---
    tasks = getTasks(config, 'chemdyg')
    if (len(tasks) == 0):
        return existing_bundles

    # --- Generate and submit chemdyg scripts ---
    for c in tasks:

        if 'ts_num_years' in c.keys():
          c['ts_num_years'] = int(c['ts_num_years'])

        # Component
         # c['component'] = getComponent(c['input_files'])
         # c['component2'] = getComponent(c['input_files2'])

        # Loop over year sets
        year_sets = getYears(c['years'])
        for s in year_sets:
            c['year1'] = s[0]
            c['year2'] = s[1]
            c['ypf'] = s[1] - s[0] + 1
            c['scriptDir'] = scriptDir
            if c['subsection'] == "ts_diags":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_ts.bash' )
            elif c['subsection'] == "o3_hole_diags":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_o3_hole_diags.bash' )
            elif c['subsection'] == "TOZ_eq_native":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_TOZ_eq_plot.bash' )
            elif c['subsection'] == "surf_o3_diags":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_surf_o3_diags.bash' )
            elif c['subsection'] == "STE_flux_native":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_STE_flux_diags.bash' )
            elif c['subsection'] == "temperature_eq_native":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_temperature_eq_plot.bash' )
            elif c['subsection'] == "summary_table_native":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_summary_table.bash' )
            elif c['subsection'] == "cmip_comparison":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_cmip_comparison.bash' )
            elif c['subsection'] == "noaa_co_comparison":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_noaa_co_comparison.bash' )
            elif c['subsection'] == "pres_lat_plots":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_pres_lat_plots.bash' )
            elif c['subsection'] == "lat_lon_plots":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_lat_lon_plots.bash' )
            elif c['subsection'] == "nox_emis_plots":
                sub = c['subsection']
                prefix = 'chemdyg_%s_%04d-%04d' % (sub,c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_nox_emis_plots.bash' )
            elif c['subsection'] == "index":
                prefix = 'chemdyg_index_%04d-%04d' % (c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_index.bash' )
            else:
                sub = c['grid']
                prefix = 'chemdyg_diags_%04d-%04d' % (c['year1'],c['year2'])
                template = templateEnv.get_template( 'chemdyg_diags.bash' )
            print(prefix)
            c['prefix'] = prefix
            scriptFile = os.path.join(scriptDir, '%s.bash' % (prefix))
            statusFile = os.path.join(scriptDir, '%s.status' % (prefix))
            settingsFile = os.path.join(scriptDir, "%s.settings" % (prefix))
            skip = checkStatus(statusFile)
            if skip:
                continue

            # Create script
            with open(scriptFile, 'w') as f:
                f.write(template.render( **c ))
            makeExecutable(scriptFile)

            with open(settingsFile, "w") as sf:
                p = pprint.PrettyPrinter(indent=2, stream=sf)
                p.pprint(c)
                p.pprint(s)

            # List of depensencies
            export = 'NONE'
            dependencies = []
            if c['subsection'] == "cmip_comparison":
                dependencies = [ os.path.join(scriptDir, 'ts_atm_monthly_180x360_aave_%04d-%04d-%04d.status' % (c['year1'],c['year2'],c['ypf'])), ]
            elif c['subsection'] == "noaa_co_comparison":
                dependencies = [ os.path.join(scriptDir, 'ts_atm_monthly_180x360_aave_%04d-%04d-%04d.status' % (c['year1'],c['year2'],c['ypf'])), ]
            elif c['subsection'] == "surf_o3_diags":
                dependencies = [ os.path.join(scriptDir, 'ts_atm_hourly_US1.0x1.0_nco_%04d-%04d-%04d.status' % (c['year1'],c['year2'],c['ypf'])), ] 
            elif c['subsection'] == "climo_diags":
                dependencies = [ os.path.join(scriptDir, 'climo_native_aave_%04d-%04d.status' % (c['year1'],c['year2'])), ]
            elif c['subsection'] == "pres_lat_plots":
                dependencies = [ os.path.join(scriptDir, 'climo_180x360_aave_%04d-%04d.status' % (c['year1'],c['year2'])), ]
            elif c['subsection'] == "lat_lon_plots":
                dependencies = [ os.path.join(scriptDir, 'climo_180x360_aave_%04d-%04d.status' % (c['year1'],c['year2'])), ]
            elif c['subsection'] == "nox_emis_plots":
                dependencies = [ os.path.join(scriptDir, 'climo_180x360_aave_%04d-%04d.status' % (c['year1'],c['year2'])), ]

            if not c["dry_run"]:
                if c["bundle"] == "":
                    # Submit job
                    submitScript(scriptFile, statusFile, export, job_ids_file, dependFiles=dependencies)
                else:
                    print("...adding to bundle '%s'" % (c["bundle"]))

    return existing_bundles
