from collections import defaultdict  # ✅ This must be here

# Standard library imports first
import time
import os
import json
from datetime import datetime

# Third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Brightway imports
import brightway2 as bw
import bw2io as bi
import bw2calc as bc
from bw2data import Database, databases as bw_databases, utils as bw2data_utils
from bw2data import calculation_setups
import bw2data as bd

class BurdenFreeAnalyzer:
    """
    Main class for analyzing burden-free activities and modifying the intervention matrix.
    This class can be modified in function of the technosphere activities to be characterised.
    """

    def __init__(self, project_name="My_BurdenFree_Project", technosphere_database_name="technosphere_database_name", ecoinvent_version = "3.11", ecoinvent_systemmodel = "cut-off"):
        self.project_name = project_name
        self.technosphere_database_name = technosphere_database_name
        self.ecoinvent_version = ecoinvent_version
        self.ecoinvent_systemmodel = ecoinvent_systemmodel
        self.all_burden_free_activities = {}
        self.setup_complete = False

    def setup_project(self):
        """Set up a new Brightway2 project and import ecoinvent if missing"""
        if self.project_name not in bd.projects:
            bd.projects.create_project(self.project_name)
       
        bd.projects.set_current(self.project_name)
        print(f"The project '{self.project_name}' is set as current.")
        # the principles are applied to the most used LCA database, ecoinvent [https://ecochain.com/blog/lci-databases-in-lca/]
        # the cutoff system model is the simplest to understand [https://support.ecoinvent.org/system-models-1]
        if self.technosphere_database_name not in bd.databases: 
            print("Importing the specified technosphere LCA database...")
            # For ecoinvent
            username = os.getenv("EI_USERNAME")
            password = os.getenv("EI_PASSWORD")

            bi.import_ecoinvent_release(self.ecoinvent_version, self.ecoinvent_systemmodel, username, password) # change the version, and system model

        self.setup_complete = True

    def save_results(self, filename="burden_free_results.json"):
        # Save burden-free activities to a JSON file.
        with open(filename, "w") as f:
            json.dump(
                {db_name: [act.as_dict() for act in acts]
                 for db_name, acts in self.all_burden_free_activities.items()},
                f,
                indent=2
            )
        
    def load_results(self, filename="burden_free_results.json"):
        """Load burden-free activities from a JSON file."""
        if not os.path.exists(filename):
            return False
        with open(filename, "r") as f:
            data = json.load(f)
        self.all_burden_free_activities = {}
        for db_name, acts in data.items():
            self.all_burden_free_activities[db_name] = []
            for act in acts:
                try:
                    # Use (database, code) as the key
                    key = (db_name, act["code"])
                    activity = bd.get_activity(key)
                    self.all_burden_free_activities[db_name].append(activity)
                except Exception as e:
                    print(f"Error loading activity {act['code']} in {db_name}: {e}")
        return True


    @staticmethod
    def is_burden_free_activity(activity):
        """Check if an activity is burden-free"""
        if len(list(activity.production())) != 1:
            return False
        if len(list(activity.technosphere())) > 0:
            return False
        if len(list(activity.biosphere())) > 0:
            return False
        return True

    def analyze_database(self, db_name):
        """Analyze one database for burden-free activities"""
        db = Database(db_name)
        burden_free_activities = [act for act in db if self.is_burden_free_activity(act)]

        print(f"\nAnalyzing database: {db_name}")
        print(f"Found {len(burden_free_activities)} burden-free activities.")

        for act in burden_free_activities[:5]:
            print(f"- {act.get('name')} ({act.get('location')}, {act.get('unit')})")
        if len(burden_free_activities) > 5:
            print("...")

        return burden_free_activities

    def analyze_all_databases(self):
        """Analyze all databases, with caching to avoid re-analysis."""
        if not self.setup_complete:
            self.setup_project()
        # Try to load cached results
        if self.load_results():
            print("✅ Loaded burden-free activities from cache.")
            return self.all_burden_free_activities
            
        # Otherwise, perform full analysis
        print("\n--- Analyzing all databases (this may take a while) ---")    
        self.all_burden_free_activities = {}
        total = 0
        start = time.time()

        print("\n--- Analyzing all databases ---")
        for db_name in bw_databases:
            activities = self.analyze_database(db_name)
            self.all_burden_free_activities[db_name] = activities
            total += len(activities)
            
        # Save results to cache
        self.save_results()
        print("\n--- Summary ---")
        for db_name, acts in self.all_burden_free_activities.items():
            print(f"Database '{db_name}': {len(acts)} burden-free activities")
        print(f"\nTotal burden-free activities: {total}")
        print(f"--- Step 1 completed in {time.time() - start:.2f} s ---")
        
        return self.all_burden_free_activities

class BiosphereFlowManager:
    """
    Manage the product flows that you want to characterised, therefore that should be duplicated as a biosphere flow 
    """

    def __init__(self, biosphere_db_name="ecoinvent-3.11-biosphere"):
        self.biosphere_db_name = biosphere_db_name
        self.biosphere_db = Database(biosphere_db_name)

    def create_or_get_biosphere_flow(self, activity, production_exchange):
        """Create or reuse a biosphere flow for a burden-free activity"""
        ref_name = production_exchange['name']
        ref_unit = production_exchange['unit']
        compartment = 'technosphere' # the flow is from the technosphere, or outside of the product system of study, the flows are from another economic environment
        subcomp = activity.get('location', 'Unknown') # same reference flow's name can come from the same activity's name, we must distinguish them furtherly 

        for flow in self.biosphere_db:
            if (flow.get('name') == ref_name and
                flow.get('unit') == ref_unit and
                flow.get('categories') == (compartment, subcomp)):
                return flow.key, False

        new_code = str(activity.key) + "_BF_" + bw2data_utils.random_string(5) # can be changed to IF (Intermediate flow) or PF (Product Flow)
        new_flow = self.biosphere_db.new_activity(code=new_code)
        new_flow.update({
            'database': self.biosphere_db_name,
            'code': new_code,
            'name': ref_name,
            'unit': ref_unit,
            'categories': (compartment, subcomp),
            'type': 'technosphere' # not a natural resource, not a emission flow, a technosphere flow 
        })
        new_flow.save()
        return new_flow.key, True

    def prepare_exchanges(self, burden_free_activities):
        """Prepare biosphere exchanges"""
        new_exchanges = []
        flows_created = 0

        for db_name, activities in burden_free_activities.items():
            print(f"\nProcessing '{db_name}':")
            for act in activities:
                prod = list(act.production())[0]
                amount, unit = prod['amount'], prod['unit']
                flow_key, created = self.create_or_get_biosphere_flow(act, prod)
                if created:
                    flows_created += 1

                new_exchanges.append((act.key, {
                    'input': flow_key,
                    'amount': amount,
                    'unit': unit,
                    'type': 'biosphere',
                    'output': act.key
                }))
                print(f"  Prepared exchange for {act.get('name')}")
        return new_exchanges, flows_created

    def add_exchanges(self, new_exchanges):
        """Add biosphere exchanges to burden-free activities"""
        for act_key, exc in new_exchanges:
            act = bw2data_utils.get_activity(act_key)
            if any(e['input'] == exc['input'] and e['type'] == exc['type'] for e in act.exchanges()):
                continue
            new_exc = act.new_exchange(**exc)
            new_exc.save()

    def process(self, burden_free_activities):
        """End-to-end biosphere flow creation"""
        print("\n--- Step 2: Creating Biosphere Flows ---")
        start = time.time()

        new_exchanges, flows_created = self.prepare_exchanges(burden_free_activities)
        print(f"Created {flows_created} new biosphere flows in '{self.biosphere_db.name}'.")

        self.add_exchanges(new_exchanges)

        print(f"--- Step 2 completed in {time.time() - start:.2f} s ---")

class CircularityCalculator:
    """
    Class to compute circularity indicators based on inventory results.
    """

    def __init__(self, excluded_flows, valuable_compartments):
        self.EXCLUDED_FLOWS = excluded_flows
        self.VALUABLE_COMPARTMENTS = valuable_compartments
        
    def should_exclude_flow(self, unit, name, categories, flow_type):
        """
        Determine if a flow should be excluded from circularity calculations.
        """
        if isinstance(categories, str):
            try:
                categories = eval(categories)
            except:
                categories = ()
                
        # 1. Excluded emissions (dependant of the list in the main notebook        
        is_excluded_emission = False
        if flow_type == 'emission' and unit == 'kilogram':
            for excluded_flow in self.EXCLUDED_FLOWS:
                if excluded_flow in name:
                    is_excluded_emission = True
                    break
                    
        # Exclude water flows expressed in m³, in ei 3.11 cutoff, this concerns 23 elementary flows - when we add the cutoff activities, 2 additional wastewater flows may be concerned: ""
        
        is_water_flow = (
            unit.lower() == 'cubic meter' 
            and 'water' in name.lower() 
            # and flow_type in ('natural resource', 'emission') # if we want to focus on initial biosphere flows           
        )
        
        return (
            unit == 'kilo Becquerel'
            or 'volume occupied' in name.lower()
            or (
                categories and len(categories) > 1
                and categories[0] == 'natural resource'
                and categories[1] == 'biotic'
                and unit in ('megajoule', 'kilowatt hour')
            )
            or is_excluded_emission
            # or is_water_flow  # if commented, the water flows are not excluded
        )

    @staticmethod
    def _subcomp_from_categories(categories):
        """
        Extract subcompartment from categories tuple.
        """
        if not categories or not isinstance(categories, tuple):
            return "unknown"

        if len(categories) > 1:
            return categories[1]
        elif len(categories) == 1:
            return categories[0]
        else:
            return "unknown"

    @staticmethod
    def build_combined_lookup(technosphere_db_name="ecoinvent-3.11-cutoff"):
        """
        Build combined lookup for both mass and energy properties in one pass.
        """
        biosphere_db = Database("ecoinvent-3.11-biosphere")
        techno_to_bio_map = {}

        for flow in biosphere_db:
            categories = flow.get('categories', ())
            if categories and categories[0] == 'technosphere' and len(categories) > 1:
                location = categories[1]
                key = (flow.get('name', ''), location)
                techno_to_bio_map[key] = (flow.get('name', ''), categories)

        mass_lookup = {}
        energy_lookup = {}
        db = Database(technosphere_db_name)

        for act in db:
            for exc in act.exchanges():
                properties = exc.get('properties', {})
                wet_mass = properties.get('wet mass', {})
                heating_value = properties.get('heating value, net', {})

                if not (wet_mass and wet_mass.get('amount') is not None) and \
                   not (heating_value and heating_value.get('amount') is not None):
                    continue

                if exc.get('type') == 'biosphere':
                    biosphere_flow = exc.input
                    if biosphere_flow is None:
                        continue
                    categories = biosphere_flow.get('categories', ())
                    cat = CircularityCalculator._subcomp_from_categories(categories)
                    key = (biosphere_flow.get('name', ''), cat)
                else:
                    flow_name = exc['name']
                    location = act['location']
                    map_key = (flow_name, location)

                    if map_key in techno_to_bio_map:
                        bio_name, categories = techno_to_bio_map[map_key]
                        cat = CircularityCalculator._subcomp_from_categories(categories)
                        key = (bio_name, cat)
                    else:
                        key = (flow_name, "technosphere")

                if wet_mass and wet_mass.get('amount') is not None:
                    mass_lookup[key] = (wet_mass.get('amount'), wet_mass.get('unit'))

                if heating_value and heating_value.get('amount') is not None:
                    energy_lookup[key] = (heating_value.get('amount'), heating_value.get('unit'))

        return mass_lookup, energy_lookup

    def get_property_factors(self, flow, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys):
        """
        Get both kg and MJ factors for a flow in one function call.
        """
        kg_factor = None
        MJ_factor = None

        unit = flow.get('unit', '')
        name = flow.get('name', '').lower()
        categories = flow.get('categories', ())
        f_type = flow.get('type', '')
        f_key = flow.key

        if self.should_exclude_flow(unit, name, categories, f_type):
            return None, None

        if unit == 'kilogram':
            kg_factor = 1.0
        elif unit == 'megajoule':
            MJ_factor = 1.0
        elif unit == 'kilowatt hour':
            MJ_factor = 3.6

        flow_props = flow.get('properties', {})
        if not kg_factor and 'wet mass' in flow_props:
            kg_factor = flow_props['wet mass'].get('amount')
        if not MJ_factor and 'heating value, net' in flow_props:
            MJ_factor = flow_props['heating value, net'].get('amount')

        if not MJ_factor and f_key in ced_cf_dict and ced_cf_dict[f_key] not in (None, 0):
            MJ_factor = ced_cf_dict[f_key]

        if f_key in renewable_cf_keys:
            MJ_factor = None

        subc = CircularityCalculator._subcomp_from_categories(categories)
        flow_name = flow.get('name', '')

        key_formats = [
            (flow_name, subc),
            (flow_name, "technosphere"),
        ]

        if len(categories) > 1:
            key_formats.append((flow_name, categories[1]))

        for key in key_formats:
            if not kg_factor and key in mass_lookup:
                kg_factor = mass_lookup[key][0]
            if not MJ_factor and key in energy_lookup:
                MJ_factor = energy_lookup[key][0]
            if kg_factor is not None and MJ_factor is not None:
                break

        return kg_factor, MJ_factor

    @staticmethod
    def get_inventory_flows(setup_name):
        """
        Retrieve inventory flows from a calculation setup.
        """
        setup = calculation_setups[setup_name]
        functional_units = setup['inv']
        all_flows = []

        for fu_dict in functional_units:
            lca = bc.LCA(fu_dict) #bw. for brightway2
            lca.lci()
            inventory = lca.inventory.sum(axis=1)
            if hasattr(inventory, 'A1'):
                inventory = inventory.A1
            reverse_dict = {v: k for k, v in lca.biosphere_dict.items()}
            non_zero_indices = np.where(np.abs(inventory) > 1e-20)[0]

            for idx in non_zero_indices:
                amount = inventory[idx]
                flow_key = reverse_dict.get(idx)
                if flow_key:
                    try:
                        flow = bd.get_activity(flow_key) #bw. for brightway2
                        all_flows.append({
                            'Flow Key': flow_key,
                            'Flow Name': flow.get('name', 'Unknown'),
                            'Amount': amount,
                            'Unit': flow.get('unit', 'Unknown'),
                            'Categories': flow.get('categories', ()),
                            'Type': flow.get('type', 'Unknown')
                        })
                    except Exception as e:
                        print(f"Error getting flow {flow_key}: {e}")

        return pd.DataFrame(all_flows) if all_flows else pd.DataFrame()

    def compute_circularity_efficiency_variables(self, flows_df, technosphere_db_name="ecoinvent-3.11-cutoff", save_csv=True, setup_name=None):
        """
        Optimized circularity computation with functional unit reference product included in Rr.
        """
        if flows_df.empty:
            return None, None, None, None

        mass_lookup, energy_lookup = CircularityCalculator.build_combined_lookup(technosphere_db_name)

        try:
            ced_method = bd.Method(('Cumulative Energy Demand (CED)', 'energy resources: non-renewable', 'energy content (HHV)')) #bw. for brightway2
            ced_cf_dict = dict(ced_method.load())
        except Exception:
            ced_cf_dict = {}

        try:
            renewable_ced_method = bd.Method(('Cumulative Energy Demand (CED)', 'energy resources: renewable', 'energy content (HHV)')) 
            renewable_cf_keys = set(k for k, _v in renewable_ced_method.load())
        except Exception:
            renewable_cf_keys = set()

        V = defaultdict(float)
        Ri = defaultdict(float)
        Rr = defaultdict(float)
        Er = defaultdict(float)
        W = defaultdict(float)
        V_mass = 0.0
        Ri_mass = 0.0
        Rr_mass = 0.0
        Er_mass = 0.0
        W_mass = 0.0
        V_energy = 0.0
        Ri_energy = 0.0
        Rr_energy = 0.0
        Er_energy = 0.0
        W_energy = 0.0
        
        detailed_flows = []
        flow_cache = {}
        
        for idx, row in flows_df.iterrows():
            flow_key = row['Flow Key']
            if flow_key not in flow_cache:
                try:
                    flow_cache[flow_key] = bd.get_activity(flow_key) 
                except Exception:
                    flow_cache[flow_key] = None

        if setup_name:
            try:
                setup = calculation_setups[setup_name]
                functional_units = setup['inv']

                for fu_dict in functional_units:
                    for fu_key, fu_amount in fu_dict.items():
                        try:
                            fu_activity = bd.get_activity(fu_key) 
                            fu_name = fu_activity.get('name', 'Unknown')
                            fu_unit = fu_activity.get('unit', 'Unknown')
                            fu_amount = abs(fu_amount)

                            print(f"Adding functional unit to Rr: {fu_amount} {fu_unit} of {fu_name}")

                            kg_factor, MJ_factor = self.get_property_factors(
                                fu_activity, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys)

                            Rr[fu_unit] += fu_amount
                            if kg_factor:
                                Rr_mass += fu_amount * kg_factor
                            if MJ_factor:
                                Rr_energy += fu_amount * MJ_factor

                            detailed_flows.append({
                                'Flow Key': fu_key,
                                'Flow Name': f"FUNCTIONAL UNIT: {fu_name}",
                                'Amount': fu_amount,
                                'Unit': fu_unit,
                                'Categories': "functional_unit",
                                'Type': "technosphere",
                                'Category': 'Recycled Outputs (Rr) - Functional Unit',
                                'Mass_Equivalent_kg': fu_amount * kg_factor if kg_factor else None,
                                'Energy_Equivalent_MJ': fu_amount * MJ_factor if MJ_factor else None,
                                'kg_Factor': kg_factor,
                                'MJ_Factor': MJ_factor
                            })

                        except Exception as e:
                            print(f"Error processing functional unit {fu_key}: {e}")
            except Exception as e:
                print(f"Error getting functional unit information: {e}")

        for idx, row in flows_df.iterrows():
            amount = float(row['Amount'])
            unit = row['Unit']
            categories = row['Categories']
            flow_type = (row['Type'] or '').lower()
            flow_key = row['Flow Key']
            name = row['Flow Name']
            flow = flow_cache[flow_key]
            if flow is None:
                continue

            if self.should_exclude_flow(unit, name, categories, flow_type):
                continue

            kg_factor, MJ_factor = self.get_property_factors(
                flow, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys)

            category = None
            mass_equiv = abs(amount) * kg_factor if kg_factor else None
            energy_equiv = abs(amount) * MJ_factor if MJ_factor else None

            if 'resource' in flow_type or any(isinstance(cat, str) and 'resource' in cat.lower() for cat in categories if categories):
                V[unit] += abs(amount)
                if kg_factor:
                    V_mass += abs(amount) * kg_factor
                if MJ_factor:
                    V_energy += abs(amount) * MJ_factor
                category = 'Natural Resources (V)'
            elif 'technosphere' in flow_type:
                if amount > 0:
                    Ri[unit] += amount
                    if kg_factor:
                        Ri_mass += amount * kg_factor
                    if MJ_factor:
                        Ri_energy += amount * MJ_factor
                    category = 'Technosphere Inputs (Ri)'
                else:
                    Rr[unit] += abs(amount)
                    if kg_factor:
                        Rr_mass += abs(amount) * kg_factor
                    if MJ_factor:
                        Rr_energy += abs(amount) * MJ_factor
                    category = 'Recycled Outputs (Rr)'
            elif 'emission' in flow_type:
                if (name == 'Water' and isinstance(categories, (list, tuple)) and len(categories) >= 2
                    and isinstance(categories[0], str) and categories[0].lower() == 'water'
                    and any(isinstance(cat, str) and cat.lower().startswith(tuple(self.VALUABLE_COMPARTMENTS)) for cat in categories)):
                    # Rr[unit] += abs(amount)
                    Er[unit] += abs(amount)
                    if kg_factor:
                        # Rr_mass += abs(amount) * kg_factor
                        Er_mass += abs(amount) * kg_factor
                    if MJ_factor:
                        # Rr_energy += abs(amount) * MJ_factor
                        Er_energy += abs(amount) * MJ_factor
                    category = 'Cleaned Emissions (Er)'
                elif (unit == "kilogram" and flow.get("name") in self.EXCLUDED_FLOWS):
                    category = 'Excluded Emission'
                    continue
                else:
                    W[unit] += abs(amount)
                    if kg_factor:
                        W_mass += abs(amount) * kg_factor
                    if MJ_factor:
                        W_energy += abs(amount) * MJ_factor
                    category = 'Waste Emissions (W)'

            detailed_flows.append({
                'Flow Key': flow_key,
                'Flow Name': name,
                'Amount': amount,
                'Unit': unit,
                'Categories': str(categories),
                'Type': flow_type,
                'Category': category,
                'Mass_Equivalent_kg': mass_equiv,
                'Energy_Equivalent_MJ': energy_equiv,
                'kg_Factor': kg_factor,
                'MJ_Factor': MJ_factor
            })

        detailed_flows_df = pd.DataFrame(detailed_flows)

        if save_csv and not detailed_flows_df.empty:
            # 1. Define the folder name
            output_folder = "results_csv"
            # 2. Create the folder if it doesn't exist
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            # 3. Combine folder and filename
            csv_filename = os.path.join(output_folder, "circularity_inventory.csv")
            detailed_flows_df.to_csv(csv_filename, index=False)
            # 4. Save
            detailed_flows_df.to_csv(csv_filename, index=False)
            print(f"Detailed inventory saved to: {os.path.abspath(csv_filename)}")

        V_df = pd.DataFrame.from_dict(V, orient='index', columns=['Natural Resources'])
        Ri_df = pd.DataFrame.from_dict(Ri, orient='index', columns=['Technosphere Inputs'])
        Rr_df = pd.DataFrame.from_dict(Rr, orient='index', columns=['Recycled Outputs'])
        Er_df = pd.DataFrame.from_dict(Er, orient='index', columns=['Cleaned Emissions'])
        W_df = pd.DataFrame.from_dict(W, orient='index', columns=['Waste Emissions'])

        combined_df = pd.concat([V_df, Ri_df, Rr_df, Er_df, W_df], axis=1).fillna(0)
        combined_df['Total Input (V+Ri)'] = combined_df['Natural Resources'] + combined_df['Technosphere Inputs']
        combined_df['Total functional (Er+Rr)'] = combined_df['Recycled Outputs'] + combined_df['Cleaned Emissions'] 
        combined_df['inefficiency (eta-)'] = combined_df['Waste Emissions'] / combined_df['Total Input (V+Ri)']
        combined_df['efficiency (eta+)'] = combined_df['Total functional (Er+Rr)'] / combined_df['Total Input (V+Ri)']

        if V_mass or Ri_mass or Rr_mass or Er_mass or W_mass:
            combined_df.loc['kilogram', ['Natural Resources', 'Technosphere Inputs', 'Recycled Outputs', 'Cleaned Emissions', 'Waste Emissions']] = [
                V_mass, Ri_mass, Rr_mass, Er_mass, W_mass
            ]
            combined_df.loc['kilogram', 'Total Input (V+Ri)'] = V_mass + Ri_mass
            combined_df.loc['kilogram', 'Total functional (Er+Rr)'] = Er_mass + Rr_mass
            combined_df.loc['kilogram', 'inefficiency (eta-)'] = (W_mass / (V_mass + Ri_mass)) if (V_mass + Ri_mass) else 0
            combined_df.loc['kilogram', 'efficiency (eta+)'] = ((Er_mass + Rr_mass) / (V_mass + Ri_mass)) if (V_mass + Ri_mass) else 0

        if V_energy or Ri_energy or Rr_energy or Er_mass or W_energy:
            combined_df.loc['megajoule', ['Natural Resources', 'Technosphere Inputs', 'Recycled Outputs', 'Cleaned Emissions', 'Waste Emissions']] = [
                V_energy, Ri_energy, Rr_energy, Er_energy, W_energy
            ]
            combined_df.loc['megajoule', 'Total Input (V+Ri)'] = V_energy + Ri_energy
            combined_df.loc['megajoule', 'Total functional (Er+Rr)'] = Er_energy + Rr_energy
            combined_df.loc['megajoule', 'inefficiency (eta-)'] = (W_energy / (V_energy + Ri_energy)) if (V_energy + Ri_energy) else 0
            combined_df.loc['megajoule', 'efficiency (eta+)'] = ((Er_energy + Rr_energy) / (V_energy + Ri_energy)) if (V_energy + Ri_energy) else 0

        return combined_df, combined_df['inefficiency (eta-)'], combined_df['efficiency (eta+)'], detailed_flows_df


    def compute_circularity_EMF_indicators(self, flows_df, technosphere_db_name="ecoinvent-3.11-cutoff", setup_name=None): # based on the Ellen Mac Arthur Method, simple version available in https://www.ellenmacarthurfoundation.org/material-circularity-indicator
        """
        Compute modified circularity indicators with formulas:
        LFI = (W + V) / (2*(Ri + V))
        CFI = (Rr + Er + Ri) / (2*(Ri + V))
        """
        if flows_df.empty:
            return None, None, None

        # Get the standard circularity results first
        results_df, LFI, CFI, detailed_flows_df = self.compute_circularity_efficiency_variables(flows_df, technosphere_db_name, save_csv=False, setup_name=setup_name)
        
        if results_df is None:
            return None, None, None

        # Calculate modified indicators for each unit
        modified_results = {}
        
        for unit in results_df.index:
            V = results_df.loc[unit, 'Natural Resources']
            Ri = results_df.loc[unit, 'Technosphere Inputs']
            Rr = results_df.loc[unit, 'Recycled Outputs']
            Er = results_df.loc[unit, 'Cleaned Emissions']
            W = results_df.loc[unit, 'Waste Emissions']

            # Calculate modified LFI and CFI
            denominator = 2 * (Ri + V)
            
            if denominator != 0:
                modified_LFI = (W + V) / denominator
                modified_CFI = (Rr + Er + Ri) / denominator
            else:
                modified_LFI = 0
                modified_CFI = 0

            modified_results[unit] = {
                'simplified_LFI': modified_LFI, # linear / total use
                'CFI': modified_CFI, # circular / total use
                'Natural_Resources_V': V,
                'Technosphere_Inputs_Ri': Ri,
                'Recycled_Outputs_Rr': Rr,
                'Cleaned_Emissions_Er': Er,
                'Waste_Emissions_W': W
            }

        # Create modified results DataFrame
        modified_df = pd.DataFrame(modified_results).T
        
        return modified_df, results_df, detailed_flows_df


class MultiLCACalculator:
    """
    Enhanced class to create and run MultiLCA calculations with circularity indicators.
    Supports multiple activities, custom amounts, and manual setup creation.
    """

    def __init__(self):
        self.location_lookup = {}

    def get_circularity_methods(self):
        """Get all circularity indicator methods."""
        return [m for m in bd.methods if m and len(m) > 0 and #bw. for brightway2
                str(m[0]).lower().startswith('circularity')]

    def precompute_locations(self, database_name="ecoinvent-3.11-cutoff"):
        """Precompute locations for all activities in the database."""
        db = bd.Database(database_name) #bw. for brightway2
        self.location_lookup = {
            activity.key: activity.get('location', 'Unknown')
            for activity in db
        }
        return self.location_lookup

    def find_activities(self, search_terms, location_codes=None, database_name="ecoinvent-3.11-cutoff"):
        """
        Find activities matching search terms and location codes.
        search_terms can be a list of multiple reference products to search for.
        """
        if not self.location_lookup:
            self.precompute_locations(database_name)

        db = bd.Database(database_name) 
        location_codes = [code.lower() for code in (location_codes or [])]
        search_terms = [term.lower() for term in search_terms]

        matches = []
        for act in db:
            name = act.get('name', '').lower()
            ref_prod = act.get('reference product', '').lower()
            unit = act.get('unit', '').lower()
            location = self.location_lookup.get(act.key, "").lower()

            # Check if any search term matches name or reference product
            term_match = any(term in name or term in ref_prod for term in search_terms)
            location_match = not location_codes or any(code in location for code in location_codes)

            if term_match and location_match:
                matches.append({
                    'key': act.key,
                    'name': act.get('name', ''),
                    'reference_product': act.get('reference product', ''),
                    'unit': act.get('unit', ''),
                    'location': self.location_lookup.get(act.key, 'Unknown')
                })

        if not matches:
            print("No activities found. Here are some similar activities:")
            similar = []
            for act in db:
                name = act.get('name', '').lower()
                ref_prod = act.get('reference product', '').lower()
                if any(term in name or term in ref_prod for term in search_terms):
                    similar.append({
                        'name': act.get('name', ''),
                        'reference_product': act.get('reference product', ''),
                        'unit': act.get('unit', ''),
                        'location': self.location_lookup.get(act.key, 'Unknown')
                    })
                    if len(similar) >= 10:
                        break

            if similar:
                df = pd.DataFrame(similar)
                print(df.to_string(index=False))
            raise ValueError("No matching activities found")

        return pd.DataFrame(matches)

    def select_activities(self):
        """
        Interactive activity selection with support for multiple activities and custom amounts.
        Returns a list of (activity_key, amount) tuples.
        """
        print("\n=== Activity Selection ===")

        # Get search parameters with better instructions
        print("Enter reference products to search for (comma-separated):")
        print("Example: electricity, high voltage, steel, concrete")
        search_terms = input("> ").strip().split(',')
        search_terms = [term.strip() for term in search_terms if term.strip()]

        location_input = input("Enter location codes (comma-separated, or leave empty): ").strip()
        location_codes = [code.strip().lower() for code in location_input.split(',')] if location_input else None

        # Find matching activities
        df = self.find_activities(search_terms, location_codes)

        print(f"\nFound {len(df)} matching activities:")
        print(df.to_string(index=True))

        # Let user select activities
        selected_activities = []
        while True:
            try:
                selection = input("\nEnter activity indices to add (comma-separated), or 'done' to finish: ")
                if selection.lower() == 'done':
                    if not selected_activities:
                        print("You must select at least one activity.")
                        continue
                    break

                indices = [int(idx.strip()) for idx in selection.split(',')]
                for idx in indices:
                    if 0 <= idx < len(df):
                        activity = df.iloc[idx]
                        print(f"\nSelected: {activity['reference_product']} ({activity['unit']})")

                        # Get amount for this activity
                        while True:
                            try:
                                amount = float(input(f"Enter amount for this activity (unit: {activity['unit']}): "))
                                if amount <= 0:
                                    print("Amount must be positive.")
                                    continue
                                selected_activities.append((activity['key'], amount))
                                break
                            except ValueError:
                                print("Please enter a valid number.")
                    else:
                        print(f"Index {idx} is out of range (0-{len(df)-1}).")
            except ValueError:
                print("Please enter valid indices or 'done'.")

        return selected_activities

    def create_calculation_setup_interactive(self, setup_name):
        """
        Create a calculation setup interactively with multiple activities.
        """
        print(f"\n=== Creating Calculation Setup: {setup_name} ===")

        # Select activities
        functional_units = self.select_activities()

        # Get circularity methods
        circ_methods = self.get_circularity_methods()

        if not circ_methods:
            raise ValueError("No circularity methods found. Available methods should start with 'Circularity'")

        print(f"\nFound {len(circ_methods)} circularity methods:")
        for i, method in enumerate(circ_methods[:5]):  # Show first 5 as examples
            print(f"{i+1}. {' | '.join(str(m) for m in method)}")
        if len(circ_methods) > 5:
            print(f"... and {len(circ_methods)-5} more")

        # Create calculation setup
        calculation_setups[setup_name] = {
            "inv": [dict([fu]) for fu in functional_units],
            "ia": circ_methods,
            "description": f"MultiLCA with {len(functional_units)} functional units and {len(circ_methods)} circularity methods"
        }

        print(f"\nCalculation setup '{setup_name}' created successfully!")
        print(f"- Functional units: {len(functional_units)}")
        print(f"- Methods: {len(circ_methods)} circularity indicators")
        for i, (key, amount) in enumerate(functional_units):
            act = bd.get_activity(key)
            print(f"  {i+1}. {act.get('reference product')} ({act.get('unit')}): {amount}")

        return calculation_setups[setup_name]

    def create_calculation_setup_manual(self, setup_name):
        """
        Create a calculation setup manually by entering JSON-like data.
        """
        print(f"\n=== Manual Calculation Setup Creation: {setup_name} ===")
        print("Enter your functional units in the format: activity_key:amount (one per line)")
        print("Example: ('database', 'key'):1.5")
        print("Enter 'done' when finished.")

        functional_units = []
        while True:
            fu_input = input(f"Functional unit {len(functional_units)+1}: ").strip()
            if fu_input.lower() == 'done':
                if not functional_units:
                    print("You must enter at least one functional unit.")
                    continue
                break
            try:
                # Parse the input
                key_str, amount_str = fu_input.split(':')
                key = eval(key_str.strip())  # Convert string to tuple
                amount = float(amount_str.strip())

                # Validate the activity exists
                bd.get_activity(key)
                functional_units.append((key, amount))
            except Exception as e:
                print(f"Invalid input: {str(e)}. Please try again.")

        # Get circularity methods
        circ_methods = self.get_circularity_methods()
        if not circ_methods:
            raise ValueError("No circularity methods found")

        # Create calculation setup
        calculation_setups[setup_name] = {
            "inv": [dict([fu]) for fu in functional_units],
            "ia": circ_methods,
            "description": f"Manual MultiLCA setup with {len(functional_units)} functional units"
        }

        print(f"\nCalculation setup '{setup_name}' created successfully!")
        print(f"- Functional units: {len(functional_units)}")
        print(f"- Methods: {len(circ_methods)} circularity indicators")
        for i, (key, amount) in enumerate(functional_units):
            act = bd.get_activity(key)
            print(f"  {i+1}. {act.get('reference product')} ({act.get('unit')}): {amount}")

        return calculation_setups[setup_name]

    def run_multilca(self, setup_name):
        """
        Run MultiLCA for a given setup and return results.
        """
        if setup_name not in calculation_setups:
            raise ValueError(f"Calculation setup '{setup_name}' not found. Available setups: {list(calculation_setups.keys())}")

        # Create MultiLCA object
        mLCA = bc.MultiLCA(setup_name) #bw. for brightway2

        # Create results DataFrame
        calculation_setup = calculation_setups[setup_name]

        columns = []
        for fu_item in calculation_setup['inv']:
            for activity_key, amount in fu_item.items():
                activity = bd.get_activity(activity_key)
                columns.append(f"{activity['name']} ({activity['location']}) - {amount} {activity.get('unit', 'unit')}")

        method_names = []
        for method in calculation_setup['ia']:
            if isinstance(method, tuple):
                method_names.append(" | ".join(str(m) for m in method))
            else:
                method_names.append(str(method))

        # Create and return the DataFrame
        mLCAdf = pd.DataFrame(
            index=method_names,
            columns=columns,
            data=mLCA.results.T
        )

        return mLCAdf

    def plot_results2(self, results_df, title=None):
        """
        Plot the MultiLCA results DataFrame grouped by LCIA categories with units.
        """
        if not isinstance(results_df, pd.DataFrame):
            raise ValueError("Input must be a DataFrame from run_multilca()")

        # Group methods by category
        material_with_water = []
        material_without_water = []
        energy = []

        # Extract units from method names
        units = {}

        for method in results_df.index:
            if isinstance(method, tuple):
                # Method is in tuple format: ("Circularity Indicator", "Material flows", "Mass - Virgin Resource Input", "kg-eq")
                if len(method) >= 4:
                    category = method[1]  # Second element is category (Material flows/Energy flows)
                    method_name = method[2]  # Third element is method name
                    unit = method[3]  # Fourth element is unit (kg-eq/MJ-eq)

                    units[method] = unit

                    if "Material flows" in category:
                        if "Without Water" in method_name:
                            material_without_water.append(method)
                        else:
                            material_with_water.append(method)
                    elif "Energy flows" in category:
                        energy.append(method)
            elif isinstance(method, str):
                # Fallback for string methods (less likely in your case)
                method_parts = method.split(" | ")
                if len(method_parts) >= 3:
                    category = method_parts[1]
                    method_name = method_parts[2]

                    # Try to extract unit from method name
                    if "kg-eq" in method_name:
                        unit = "kg-eq"
                    elif "MJ-eq" in method_name:
                        unit = "MJ-eq"
                    else:
                        unit = "unit"

                    units[method] = unit

                    if "Material flows" in category:
                        if "Without Water" in method_name:
                            material_without_water.append(method)
                        else:
                            material_with_water.append(method)
                    elif "Energy flows" in category:
                        energy.append(method)

        # Create a figure with subplots for each category
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 18))
        if title:
            fig.suptitle(title, y=1.02)
        else:
            fig.suptitle("MultiLCA Results by Category", y=1.02)

        # Plot each category
        categories = [
            ("Material Flows (with water)", material_with_water, axes[0]),
            ("Material Flows (without water)", material_without_water, axes[1]),
            ("Energy Flows", energy, axes[2])
        ]

        for category_name, methods, ax in categories:
            if methods:
                # Filter dataframe for this category
                category_df = results_df.loc[methods]

                # Plot
                category_df.plot(kind='bar', ax=ax, legend=False)

                # Formatting
                ax.set_title(category_name)

                # Get the unit for this category
                if methods:
                    sample_method = methods[0]
                    unit = units.get(sample_method, "unit")
                    ax.set_ylabel(f"Impact Score ({unit})")

                ax.grid(axis='y', alpha=0.3)

                # Rotate x-axis labels
                plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
            else:
                ax.axis('off')
                ax.set_title(f"No {category_name} methods found")

        plt.tight_layout()
        plt.show()

    def plot_results(self, results_df, title=None):
        """Plot the MultiLCA results DataFrame."""
        plt.figure(figsize=(12, 8))
        results_df.plot(kind='bar', ax=plt.gca())

        if title:
            plt.title(title)
        else:
            plt.title("MultiLCA Results")

        plt.ylabel("Impact Score")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        plt.show()

    def full_workflow(self):
        """
        Complete workflow from setup creation to results visualization.
        """
        print("\n=== MultiLCA Calculation Workflow ===")

        try:
            # Choose creation method
            print("\nSelect setup creation method:")
            print("1. Interactive setup with activity search")
            print("2. Manual setup creation")
            choice = input("Enter your choice (1 or 2): ").strip()

            setup_name = input("\nEnter a name for this calculation setup: ").strip()

            if choice == "1":
                # Interactive setup
                self.create_calculation_setup_interactive(setup_name)
            elif choice == "2":
                # Manual setup
                self.create_calculation_setup_manual(setup_name)
            else:
                print("Invalid choice. Please enter 1 or 2.")
                return None

            # Run MultiLCA
            print("\nRunning MultiLCA...")
            results = self.run_multilca(setup_name)

            # Show and plot results
            print("\nMultiLCA Results:")
            print(results)

            plot_title = f"Results for {setup_name}"
            self.plot_results2(results, title=plot_title) # or plot_results1 for rough ploting.

            return results

        except Exception as e:
            print(f"\n❌ Error: {str(e)}")
            import traceback
            traceback.print_exc()
            return None

class CircularityDatabaseAnalyzer:
    """
    Class to compute circularity indicators for all activities in the database,
    classify them by location or ISIC class, and plot the results.
    """
    def __init__(self, isic_section_map, isic_division_map, location_to_un_group, excluded_flows, valuable_compartments):
        self.EXCLUDED_FLOWS = excluded_flows
        self.VALUABLE_COMPARTMENTS = valuable_compartments
        self.ISIC_SECTION_MAP = isic_section_map
        self.ISIC_DIVISION_MAP = isic_division_map
        self.LOCATION_TO_UN_GROUP = location_to_un_group
        self.multi_lca_calculator = MultiLCACalculator() #location
        self.circularity_calculator = CircularityCalculator(excluded_flows, valuable_compartments)

    def get_isic_class(self, activity):
        """
        Retrieve ISIC class from Brightway2 activity metadata.
        """
        isic_class = activity.get('ISIC') or activity.get('ISIC code') or 'Unknown'
        if isic_class == 'Unknown':
            classifications = activity.get('classifications', [])
            for classification in classifications:
                if isinstance(classification, tuple) and classification[0] == 'ISIC rev.4 ecoinvent':
                    isic_class = classification[1]
                    break
        return isic_class

    def parse_isic(self, isic_code):
        """
        Parse ISIC class into section, division, and description.
        """
        try:
            division = isic_code[:2]
            section, description = self.ISIC_DIVISION_MAP.get(division, ("Unknown", "Unknown"))
            return {
                "division": division,
                "section": section,
                "description": description
            }
        except Exception:
            return {
                "division": "Unknown",
                "section": "Unknown",
                "description": "Unknown"
            }

    def precompute_isic_classes(self, database_name="ecoinvent-3.11-cutoff"):
        """
        Precompute ISIC classes for all activities in the database.
        """
        db = bd.Database(database_name) 
        isic_lookup = {}
        for activity in db:
            isic_class = self.get_isic_class(activity)
            isic_lookup[activity.key] = isic_class
        return isic_lookup

    def precompute_isic_classes_classified(self, database_name="ecoinvent-3.11-cutoff"):
        """
        Precompute ISIC classes for all activities in the database, including section and division.
        """
        db = bd.Database(database_name)
        isic_lookup = {}
        for activity in db:
            isic_class = None
            for c in activity.get('classifications', []):
                if c[0].startswith("ISIC"):
                    isic_class = c[1].split(":")[0]
                    break
            if isic_class:
                parsed = self.parse_isic(isic_class)
                isic_lookup[activity.key] = {
                    "isic_code": isic_class,
                    "division": parsed["division"],
                    "section": parsed["section"],
                    "description": parsed["description"]
                }
            else:
                isic_lookup[activity.key] = {
                    "isic_code": "Unknown",
                    "division": "Unknown",
                    "section": "Unknown",
                    "description": "Unknown"
                }
        return isic_lookup

    def precompute_locations_classified(self, database_name="ecoinvent-3.11-cutoff"):
        """
        Precompute locations for all activities in the database, classified by UN Regional Group.
        """
        db = bd.Database(database_name)
        location_lookup_classified = {}
        for activity in db:
            location = activity.get('location', 'Unknown')
            un_group = self.LOCATION_TO_UN_GROUP.get(location, 'Unknown')
            location_lookup_classified[activity.key] = {
                "location": location,
                "un_group": un_group
            }
        return location_lookup_classified

    def get_inventory_for_1_act(self, fu_dict):
        """
        Get inventory flows for a single functional unit.
        """
        all_flows = []

        try:
            # Create FU dict for LCA
            lca_fu_dict = {fu_dict['activity_key']: fu_dict['amount']}
            # Initialize LCA and perform Life Cycle Inventory (LCI)
            lca = bc.LCA(lca_fu_dict)
            lca.lci()
            # Sum the inventory
            inventory = lca.inventory.sum(axis=1)
            if hasattr(inventory, 'A1'):
                inventory = inventory.A1
            # Reverse dictionary to map indices to flow keys
            reverse_dict = {v: k for k, v in lca.biosphere_dict.items()}
            # Find non-zero indices in the inventory
            non_zero_indices = np.where(np.abs(inventory) > 1e-20)[0]
            # Collect flow information for non-zero indices
            for idx in non_zero_indices:
                amount = inventory[idx]
                flow_key = reverse_dict.get(idx)
                if flow_key:
                    try:
                        flow = bd.get_activity(flow_key) # flow = bd.get_node(flow_key) may be used in the future
                        # Verify we got a proper dictionary
                        if flow is None: # checking instance (CI) to be removed when optimised
                            print(f"Warning: flow {flow_key} is not found. Key type: {type(flow_key)}") # CI
                            continue # CI
                        all_flows.append({
                            'Flow Key': flow_key,
                            'Flow Name': flow.get('name', 'Unknown'),
                            'Amount': amount,
                            'Unit': flow.get('unit', 'Unknown'),
                            'Categories': flow.get('categories', ()),
                            'Type': flow.get('type', 'Unknown')
                        })
                    except Exception as e:
                        print(f"Error getting flow {flow_key}: {e}")
                        continue # Noisy skip for errors
                        
            # ADD THE REFERENCE PRODUCT TO THE INVENTORY
            activity = bd.get_activity(fu_dict['activity_key']) 
            # if activity is None:
                # print(f"Warning: activity {fu_dict['activity_key']} is not found.")
            all_flows.append({
                'Flow Key': fu_dict['activity_key'],
                'Flow Name': activity.get('reference product', 'Unknown'),
                'Amount': -fu_dict['amount'],
                'Unit': activity.get('unit', 'Unknown'),
                'Categories': activity.get('categories', ()),
                'Type': 'technosphere'
            })

        except Exception as e:
            print(f"LCA error for activity: {e}")
            import traceback
            print(traceback.format_exc())
            
        return pd.DataFrame(all_flows) if all_flows else pd.DataFrame()

    def compute_circularity_for_single_fu(self, fu_dict, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys, isic_lookup_classified, location_lookup_classified):
        """
        Compute circularity for a single functional unit.
        """
        try:
            activity = bd.get_activity(fu_dict['activity_key'])
            flows_df = self.get_inventory_for_1_act(fu_dict)
            if flows_df.empty:
                return None

            V = defaultdict(float)
            Ri = defaultdict(float)
            Rr = defaultdict(float)
            Er = defaultdict(float)
            W = defaultdict(float)
            V_mass = 0.0
            Ri_mass = 0.0
            Rr_mass = 0.0
            Er_mass = 0.0
            W_mass = 0.0
            V_energy = 0.0
            Ri_energy = 0.0
            Rr_energy = 0.0
            Er_energy = 0.0
            W_energy = 0.0

            flow_cache = {}
            
            for idx, row in flows_df.iterrows():
                flow_key = row['Flow Key']
                if flow_key not in flow_cache:
                    try:
                        flow = bd.get_activity(flow_key)
                        flow_cache[flow_key] = flow
                    except Exception:
                        flow_cache[flow_key] = None
                        # flow_cache[flow_key] = bw.get_activity(flow_key # may work directly
                        if flow is None:
                            continue

            for idx, row in flows_df.iterrows():
                amount = float(row['Amount'])
                unit = row['Unit']
                categories = row['Categories']
                flow_type = (row['Type'] or '').lower()
                flow_key = row['Flow Key']
                name = row['Flow Name']
                flow = flow_cache[flow_key]
                if flow is None:
                    continue
                # Safe access to flow properties
                flow_name = flow.get('name', 'Unknown')
                flow_unit = flow.get('unit', 'Unknown')
                flow_categories = flow.get('categories', ())
                
                if self.circularity_calculator.should_exclude_flow(unit, name, categories, flow_type):
                    continue

                kg_factor, MJ_factor = self.circularity_calculator.get_property_factors(
                    flow, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys)
                
                # Classification logic (V, Ri, etc.)
                if 'resource' in flow_type or any(isinstance(cat, str) and 'resource' in cat.lower() for cat in categories if categories):
                    V[unit] += abs(amount)
                    if kg_factor:
                        V_mass += abs(amount) * kg_factor
                    if MJ_factor:
                        V_energy += abs(amount) * MJ_factor
                elif 'technosphere' in flow_type:
                    if amount > 0:
                        Ri[unit] += amount
                        if kg_factor:
                            Ri_mass += amount * kg_factor
                        if MJ_factor:
                            Ri_energy += amount * MJ_factor
                    else:
                        Rr[unit] += abs(amount)
                        if kg_factor:
                            Rr_mass += abs(amount) * kg_factor
                        if MJ_factor:
                            Rr_energy += abs(amount) * MJ_factor
                elif 'emission' in flow_type:
                    if (name == 'Water' and isinstance(categories, (list, tuple)) and len(categories) >= 2
                        and isinstance(categories[0], str) and categories[0].lower() == 'water'
                        and any(isinstance(cat, str) and cat.lower().startswith(tuple(self.VALUABLE_COMPARTMENTS)) for cat in categories)):
                        Er[unit] += abs(amount)
                        if kg_factor:
                            Er_mass += abs(amount) * kg_factor
                        if MJ_factor:
                            Er_energy += abs(amount) * MJ_factor
                    # elif (unit == "kilogram" and flow.get("name") in self.EXCLUDED_FLOWS):
                        # continue
                    elif (unit == "kilogram" and any(excluded_flow.split(",")[0].strip() in str(flow.get("name", "")) for excluded_flow in self.EXCLUDED_FLOWS)):
                        continue
                    else:
                        W[unit] += abs(amount)
                        if kg_factor:
                            W_mass += abs(amount) * kg_factor
                        if MJ_factor:
                            W_energy += abs(amount) * MJ_factor

            isic_info = isic_lookup_classified.get(activity.key, {
                "isic_code": "Unknown",
                "division": "Unknown",
                "section": "Unknown",
                "description": "Unknown"
            })
            location_info = location_lookup_classified.get(activity.key, {})
            location_code = location_info.get("location", "Unknown")
            un_group = location_info.get("un_group", "Unknown")

            return {
                'Process Key': activity.key,
                'Process Name': activity.get('name'),
                'Location': activity.get('location', 'Unknown'),
                'Reference Product': activity.get('reference product', 'Unknown'),
                'Unit': activity.get('unit', 'Unknown'),
                'ISIC Code': isic_info["isic_code"],
                'ISIC Division': isic_info["division"],
                'ISIC Section': isic_info["section"],
                'ISIC Description': isic_info["description"],
                'LOCATION Code': location_code,
                'UN Group': un_group,
                'V_kg': V_mass,
                'Ri_kg': Ri_mass,
                'Rr_kg': Rr_mass,
                'Er_kg': Er_mass,
                'W_kg': W_mass,
                'eta-_kg': (W_mass / (V_mass + Ri_mass)) if (V_mass + Ri_mass) else 0, # efficiency negative perspective
                'eta+_kg': ((Rr_mass + Er_mass) / (V_mass + Ri_mass)) if (V_mass + Ri_mass) else 0,
                'LFI_kg': ((W_mass + V_mass) / (2*(V_mass + Ri_mass))) if (V_mass + Ri_mass) else 0, # EMF perspective
                'CFI_kg': ((Rr_mass + Er_mass + Ri_mass) / (2*(V_mass + Ri_mass))) if (V_mass + Ri_mass) else 0, # EMF positive perspective
                'V_MJ': V_energy,
                'Ri_MJ': Ri_energy,
                'Rr_MJ': Rr_energy,
                'Er_MJ': Er_energy,
                'W_MJ': W_energy,
                'eta-_MJ': (W_energy / (V_energy + Ri_energy)) if (V_energy + Ri_energy) else 0,
                'eta+_MJ': ((Rr_energy + Er_energy) / (V_energy + Ri_energy)) if (V_energy + Ri_energy) else 0,
                'LFI_MJ': ((W_energy + V_energy) / (2*(V_energy + Ri_energy))) if (V_energy + Ri_energy) else 0, # EMF perspective
                'CFI_MJ': ((Rr_energy + Er_energy + Ri_energy) / (2*(V_energy + Ri_energy))) if (V_energy + Ri_energy) else 0, # EMF positive perspective
            }

        except Exception as e:
            print(f"Error processing {activity.get('name') if 'activity' in locals() else 'unknown'}: {e}")
            print(f"Error processing {activity.get('name', 'unknown')}: {e}")
            import traceback
            print(traceback.format_exc())
            return None

    def compute_circularity_for_all_processes_sequential(self, database_name="ecoinvent-3.11-cutoff", sample_size=None, save_interval=500):
        """
        Compute circularity metrics for all processes sequentially.
        """
        db = bd.Database(database_name) #bw. for brightway2
        # functional_units = [] # here in the functional file
        activities = list(db)

        if not activities:
            print("No activities found in the database.")
            return pd.DataFrame()

        if sample_size and sample_size < len(activities):
            activities = pd.Series(activities).sample(sample_size, random_state=42).tolist()

        functional_units = []
        for activity in activities:
            production_exchanges = [exc for exc in activity.production() if exc.get('type') == 'production']
            if production_exchanges:
                production_exc = production_exchanges[0]
                functional_units.append({
                    'activity_key': activity.key,
                    'amount': production_exc['amount'],
                    'activity': activity
                })

        if not functional_units:
            print("No functional units prepared. Check if activities have production exchanges.")
            return pd.DataFrame()

        print(f"Prepared {len(functional_units)} functional units")

        try:
            ced_method = bd.Method(('Cumulative Energy Demand (CED)', 'energy resources: non-renewable', 'energy content (HHV)'))
            ced_cf_dict = dict(ced_method.load())
        except Exception:
            ced_cf_dict = {}

        try:
            renewable_ced_method = bd.Method(('Cumulative Energy Demand (CED)', 'energy resources: renewable', 'energy content (HHV)'))
            renewable_cf_keys = set(k for k, _v in renewable_ced_method.load())
        except Exception:
            renewable_cf_keys = set()

        mass_lookup, energy_lookup = CircularityCalculator.build_combined_lookup(database_name)
        isic_lookup_classified = self.precompute_isic_classes_classified(database_name)
        location_lookup_classified = self.precompute_locations_classified(database_name) # self.multi_lca_calculator.precompute_locations(database_name) # function already in mlca class # need to classified it

        results = []
        start_time = time.time() 

        print(f"Starting sequential processing of {len(functional_units)} processes...")

        for i, fu_dict in enumerate(functional_units):
            try:
                if i % 10 == 0 or i == len(functional_units) - 1:
                    elapsed = time.time() - start_time
                    rate = (i + 1) / elapsed if elapsed > 0 else 0
                    remaining = (len(functional_units) - i - 1) / rate if rate > 0 else float('inf')
                    print(f"Processing {i+1}/{len(functional_units)} ({((i+1)/len(functional_units))*100:.1f}%) - "
                          f"Rate: {rate:.1f} processes/min - "
                          f"ETA: {remaining/60:.1f} minutes")

                result = self.compute_circularity_for_single_fu(fu_dict, mass_lookup, energy_lookup, ced_cf_dict, renewable_cf_keys, isic_lookup_classified, location_lookup_classified)
                if result:
                    results.append(result)

                if len(results) % save_interval == 0:
                    interim_df = pd.DataFrame(results)
                    # 1. Define and create the folder
                    output_folder = "results_csv"
                    if not os.path.exists(output_folder):
                        os.makedirs(output_folder)
                    # 2. Construct the path for the interim file
                    interim_filename = os.path.join(output_folder, f'circularity_interim_{len(results)}.csv')

                    # 3. Save
                    interim_df.to_csv(interim_filename, index=False)
                    print(f"Saved interim results to {interim_filename}: {len(results)} processes")

            except Exception as e:
                print(f"Error processing activity {i+1}: {e}")
                continue

        total_time = time.time() - start_time
        print(f"Sequential processing completed in {total_time/60:.1f} minutes")

        if len(functional_units) > 0:
            success_rate = (len(results) / len(functional_units)) * 100
            print(f"Success rate: {len(results)}/{len(functional_units)} ({success_rate:.1f}%)")
        else:
            print("No functional units were processed.")

        return pd.DataFrame(results)

    def plot_hybrid_scale_scatter(self, results_df, save_path="basic_plots", group_by="section"):
        """
        Create scatter plots of LFI, CFI, eta+, and eta- in kg vs MJ with a hybrid scale.
        """
        import os
        import seaborn as sns
        import matplotlib.ticker as ticker

        if results_df.empty:
            print("No data to plot.")
            return

        if group_by.lower() == "section":
            col = "ISIC Section"
        elif group_by.lower() == "division":
            col = "ISIC Division"
        elif group_by.lower() == "description":
            col = "ISIC Description"
        elif group_by.lower() == "location":
            col = "LOCATION Code"
        elif group_by.lower() == "un_group":
            col = "UN Group"
        else:
            raise ValueError("group_by must be 'section', 'division', or 'description'")

        os.makedirs(save_path, exist_ok=True)
        unique_groups = results_df[col].dropna().unique()

        # Create a FIXED color palette that will be consistent across all plots
        # Using the same categories as in your image for consistency
        fixed_colors = {
            'A': '#1f77b4',  # Agriculture, forestry and fishing - Blue
            'B': '#ff7f0e',  # Mining and quarrying - Orange
            'C': '#2ca02c',  # Manufacturing - Green
            'D': '#d62728',  # Electricity, gas, steam - Red
            'E': '#9467bd',  # Water supply; sewerage - Purple
            'F': '#8c564b',  # Construction - Brown
            'G': '#e377c2',  # Wholesale and retail trade - Pink
            'I': '#7f7f7f',  # Accommodation and food service - Gray
            'J': '#bcbd22',  # Information and communication - Yellow-green
            'M': '#17becf',  # Professional, scientific - Cyan
            'N': '#ff9896',  # Administrative and support - Light red
            'S': '#c5b0d5'   # Other service activities - Light purple
        }
        color_palette = sns.color_palette("tab20", len(unique_groups))

        # For any groups not in the fixed mapping, generate consistent colors
        # color_map = dict(zip(unique_groups, color_palette))
        color_map = {}

        for i, grp in enumerate(unique_groups):
            if grp in fixed_colors:
                color_map[grp] = fixed_colors[grp]
            else:
                color_map[grp] = color_palette[i]

        class LinearLogTransform:
            input_dims = output_dims = 1
            def transform(self, a):
                a = np.clip(a, 1e-10, None)
                return np.where(a <= 1, a, np.log10(a) + 1)
            def inverted(self):
                return InvertedLinearLogTransform()

        class InvertedLinearLogTransform:
            input_dims = output_dims = 1
            def transform(self, a):
                return np.where(a <= 1, a, 10**(a - 1))
            def inverted(self):
                return LinearLogTransform()

        linear_log_transform = LinearLogTransform()
        inverted_linear_log_transform = InvertedLinearLogTransform()

        # LFI scatter plot
        fig1, ax1 = plt.subplots(figsize=(10, 8))
        ax1.set_xscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        ax1.set_yscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        for grp in unique_groups:
            subset = results_df[results_df[col] == grp]
            if col == "ISIC Section":
                label = f"{grp} — {self.ISIC_SECTION_MAP.get(grp, 'Unknown')}"
            else:
                label = f"{grp}"
            ax1.scatter(
                subset['LFI_kg'],
                subset['LFI_MJ'],
                color=color_map[grp],
                label=label,
                alpha=0.6,
                s=10  # Reduced from default (try values between 10-20)
            )
        ax1.xaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
        ax1.yaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
        ax1.plot([0, 1, 2, 3], [0, 1, 2, 3], 'r--', alpha=0.5, label='y = x (Equality Line)')
        ax1.set_xlabel('LFI (kg-based)')
        ax1.set_ylabel('LFI (MJ-based)')
        ax1.set_title(f'Linear Flow Index: kg vs MJ (Hybrid Scale) by {col}')
        ax1.grid(True, alpha=0.3)
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left') # , fontsize=8)
        fig1.tight_layout()
        # fig1.subplots_adjust(right=0.8)  # explicitly reserve space for the legend
        fig1.savefig(f"{save_path}/LFI_hybrid_scale_by_{group_by}.png", dpi=300, bbox_inches='tight') #400 dpi mat be suited

        # CFI scatter plot
        fig2, ax2 = plt.subplots(figsize=(10, 8))
        ax2.set_xscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        ax2.set_yscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        for grp in unique_groups:
            subset = results_df[results_df[col] == grp]
            if col == "ISIC Section":
                label = f"{grp} — {self.ISIC_SECTION_MAP.get(grp, 'Unknown')}"
            else:
                label = f"{grp}"
            ax2.scatter(
                subset['CFI_kg'],
                subset['CFI_MJ'],
                color=color_map[grp],
                label=label,
                alpha=0.6,
                s=10  
            )
        ax2.xaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
        ax2.yaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
        ax2.plot([0, 1, 2, 3], [0, 1, 2, 3], 'r--', alpha=0.5, label='y = x (Equality Line)')
        ax2.set_xlabel('CFI (kg-based)')
        ax2.set_ylabel('CFI (MJ-based)')
        ax2.set_title(f'Circular Flow Index: kg vs MJ (Hybrid Scale) by {col}')
        ax2.grid(True, alpha=0.3)
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig2.tight_layout()
        # fig2.subplots_adjust(right=0.8)  # explicitly reserve space for the legend
        fig2.savefig(f"{save_path}/CFI_hybrid_scale_by_{group_by}.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)

        # eta- scatter plot
        fig3, ax3 = plt.subplots(figsize=(10, 8))
        ax3.set_xscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        ax3.set_yscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        for grp in unique_groups:
            subset = results_df[results_df[col] == grp]
            if col == "ISIC Section":
                label = f"{grp} — {self.ISIC_SECTION_MAP.get(grp, 'Unknown')}"
            else:
                label = f"{grp}"
            ax3.scatter(
                subset['eta-_kg'],
                subset['eta-_MJ'],
                color=color_map[grp],
                label=label,
                alpha=0.6,
                s=10  
            )
        ax3.xaxis.set_major_locator(ticker.FixedLocator([0, 0.5, 1]))
        ax3.yaxis.set_major_locator(ticker.FixedLocator([0, 0.5, 1]))
        ax3.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='y = x (Equality Line)')
        ax3.set_xlabel('eta- (kg-based)')
        ax3.set_ylabel('eta- (MJ-based)')
        ax3.set_title(f'Efficiency Negative: kg vs MJ (Hybrid Scale) by {col}')
        ax3.grid(True, alpha=0.3)
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig3.tight_layout()
        # fig3.subplots_adjust(right=0.8)  # explicitly reserve space for the legend
        fig3.savefig(f"{save_path}/eta-_hybrid_scale_by_{group_by}.png", dpi=300, bbox_inches='tight')
        plt.close(fig3)

        # eta+ scatter plot
        fig4, ax4 = plt.subplots(figsize=(10, 8))
        ax4.set_xscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        ax4.set_yscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
        for grp in unique_groups:
            subset = results_df[results_df[col] == grp]
            if col == "ISIC Section":
                label = f"{grp} — {self.ISIC_SECTION_MAP.get(grp, 'Unknown')}"
            else:
                label = f"{grp}"
            ax4.scatter(
                subset['eta+_kg'],
                subset['eta+_MJ'],
                color=color_map[grp],
                label=label,
                alpha=0.6,
                s=10  
            )
        ax4.xaxis.set_major_locator(ticker.FixedLocator([0, 0.5, 1]))
        ax4.yaxis.set_major_locator(ticker.FixedLocator([0, 0.5, 1]))
        ax4.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='y = x (Equality Line)')
        ax4.set_xlabel('eta+ (kg-based)')
        ax4.set_ylabel('eta+ (MJ-based)')
        ax4.set_title(f'Efficiency Positive: kg vs MJ (Hybrid Scale) by {col}')
        ax4.grid(True, alpha=0.3)
        ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig4.tight_layout()
        # fig4.subplots_adjust(right=0.8)  # explicitly reserve space for the legend
        fig4.savefig(f"{save_path}/eta+_hybrid_scale_by_{group_by}.png", dpi=300, bbox_inches='tight')
        plt.close(fig4)

    def plot_cross_indicator_scatter(self, results_df, save_path="basic_plots", group_by="section"):
        """
        Create hybrid-scale scatter plots comparing indicators (LFI vs CFI, eta+ vs eta-) 
        for both kg and MJ bases, grouped by ISIC or location.
        """
        import os
        import seaborn as sns
        import matplotlib.ticker as ticker
        import matplotlib.pyplot as plt
        import numpy as np

        if results_df.empty:
            print("No data to plot.")
            return

        # Map grouping column
        group_cols = {
            "section": "ISIC Section",
            "division": "ISIC Division",
            "description": "ISIC Description",
            "location": "LOCATION Code",
            "un_group": "UN Group"
        }
        col = group_cols.get(group_by.lower())
        if not col:
            raise ValueError("group_by must be 'section', 'division', 'description', 'location', or 'un_group'")

        os.makedirs(save_path, exist_ok=True)
        unique_groups = results_df[col].dropna().unique()
        color_palette = sns.color_palette("tab20", len(unique_groups))
        color_map = dict(zip(unique_groups, color_palette))

        # Hybrid scale transformation
        class LinearLogTransform:
            input_dims = output_dims = 1
            def transform(self, a):
                a = np.clip(a, 1e-10, None)
                return np.where(a <= 1, a, np.log10(a) + 1)
            def inverted(self):
                return InvertedLinearLogTransform()

        class InvertedLinearLogTransform:
            input_dims = output_dims = 1
            def transform(self, a):
                return np.where(a <= 1, a, 10**(a - 1))
            def inverted(self):
                return LinearLogTransform()

        linear_log_transform = LinearLogTransform()
        inverted_linear_log_transform = InvertedLinearLogTransform()

        # Define pairs to plot
        pairs = [
            ("LFI_kg", "CFI_kg", "Linear vs Circular Flow Index (kg)", "LFI (kg)", "CFI (kg)"),
            ("LFI_MJ", "CFI_MJ", "Linear vs Circular Flow Index (MJ)", "LFI (MJ)", "CFI (MJ)"),
            ("eta+_kg", "eta-_kg", "Efficiency Positive vs Negative (kg)", "η⁺ (kg)", "η⁻ (kg)"),
            ("eta+_MJ", "eta-_MJ", "Efficiency Positive vs Negative (MJ)", "η⁺ (MJ)", "η⁻ (MJ)")
        ]

        for x_var, y_var, title, x_label, y_label in pairs:
            fig, ax = plt.subplots(figsize=(10, 8))
            ax.set_xscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))
            ax.set_yscale('function', functions=(linear_log_transform.transform, inverted_linear_log_transform.transform))

            for grp in unique_groups:
                subset = results_df[results_df[col] == grp]
                label = f"{grp} — {self.ISIC_SECTION_MAP.get(grp, 'Unknown')}" if col == "ISIC Section" else f"{grp}"
                ax.scatter(subset[x_var], subset[y_var], color=color_map[grp], label=label, alpha=0.6, s=10)

            # Define plot limits
            lim = [0, 1]
            # ax.set_xlim(lim)
            # ax.set_ylim(lim)
            # Draw the balanced line: x + y = 1 → from (0,1) to (1,0)
            ax.plot([0, 1], [1, 0], 'r--', alpha=0.5, label='x + y = 1 (balanced)')
            ax.xaxis.set_major_locator(ticker.FixedLocator(lim))
            ax.yaxis.set_major_locator(ticker.FixedLocator(lim))
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(f"{title} (Hybrid Scale) by {col}")
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            fig.tight_layout()
            # fig.subplots_adjust(right=0.8)  # explicitly reserve space for the legend
            fig.savefig(f"{save_path}/{x_var}_vs_{y_var}_hybrid_by_{group_by}.png", dpi=300, bbox_inches='tight')
            plt.close(fig)
            
class ProgressTracker:
    """Track progress across runs"""

    def __init__(self, project_name):
        self.project_name = project_name
        self.progress_file = f"progress_{project_name}.json"
        self.progress_data = self.load_progress()

    def load_progress(self):
        if os.path.exists(self.progress_file):
            try:
                with open(self.progress_file, "r") as f:
                    return json.load(f)
            except:
                return {}
        return {}

    def save_progress(self):
        with open(self.progress_file, "w") as f:
            json.dump(self.progress_data, f, indent=2)

    def update_step(self, step, status, data=None, execution_time=None):
        self.progress_data[step] = {
            "status": status,
            "timestamp": datetime.now().isoformat(),
            "execution_time": execution_time,
            "data_summary": data
        }
        self.save_progress()

    def is_step_completed(self, step):
        return step in self.progress_data and self.progress_data[step]["status"] == "completed"

    def print_summary(self):
        print("\n=== Progress Summary ===")
        for step, data in self.progress_data.items():
            print(f"{step}: {data.get('status')} at {data.get('timestamp')} (time: {data.get('execution_time')})")
