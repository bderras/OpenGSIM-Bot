# Fichier : run_gmm_calculation_v8.py
# Basé sur la v7, avec une correction d'unité ciblée pour ConvertitoEtAl2012.

import sys
import os
import json
import argparse
import importlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import inspect

sys.path.insert(0, "/app/oq-engine-source")

from openquake.hazardlib import const, imt
from openquake.hazardlib.contexts import SitesContext, RuptureContext, DistancesContext
from openquake.hazardlib.gsim.base import GMPE, CoeffsTable

def get_fas_frequencies(gmm_class):
    """Lit dynamiquement les fréquences supportées par un GMM (pour FAS, EAS, DRVT)."""
    frequencies = []
    try:
        for attr_name in dir(gmm_class):
            try:
                attr = getattr(gmm_class, attr_name)
                if isinstance(attr, CoeffsTable):
                    for key in attr._coeffs.keys():
                        if hasattr(key, 'frequency'):
                            frequencies.append(float(key.frequency))
            except:
                continue
        if not frequencies:
            frequencies = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 25.0, 50.0, 100.0]
    except:
        frequencies = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 25.0, 50.0, 100.0]
    return sorted(list(set(frequencies)))

def get_supported_imts(gmm_class):
    """Méthode intelligente pour détecter automatiquement les IMTs supportés par un modèle GMM."""
    supported_imts = []
    if hasattr(gmm_class, 'DEFINED_FOR_INTENSITY_MEASURE_TYPES'):
        defined_imts = gmm_class.DEFINED_FOR_INTENSITY_MEASURE_TYPES
        for imt_type in defined_imts:
            imt_name = imt_type.__name__
            if imt_name == "SA":
                periods = get_sa_periods(gmm_class)
                for period in periods: supported_imts.append(("SA", period))
            elif imt_name == "AvgSA":
                periods = get_sa_periods(gmm_class)
                for period in periods:
                    if period > 0: supported_imts.append(("AvgSA", period))
            elif imt_name in ["EAS", "FAS", "DRVT"]:
                frequencies = get_fas_frequencies(gmm_class)
                for freq in frequencies:
                    supported_imts.append((imt_name, freq))
            elif imt_name == "SDi":
                periods, strength_ratios = get_sdi_parameters(gmm_class)
                for period in periods:
                    for ratio in strength_ratios:
                        supported_imts.append(("SDi", {"period": period, "strength_ratio": ratio}))
            else:
                supported_imts.append((imt_name, None))
    if not supported_imts:
        base_imts = ["PGA", "PGV", "PGD", "IA", "CAV", "CAV5", "MMI", "AI", "ASI", "DSI", "SI", "VSI", "LSD"]
        duration_imts = ["RSD575", "RSD595", "RSD2080", "D5_75", "D5_95"]
        for imt_name in base_imts:
            if hasattr(imt, imt_name): supported_imts.append((imt_name, None))
        for imt_name in duration_imts:
            if hasattr(imt, imt_name): supported_imts.append((imt_name, None))
        sa_periods = get_sa_periods(gmm_class)
        for period in sa_periods: supported_imts.append(("SA", period))
    return supported_imts

def get_sdi_parameters(gmm_class):
    periods, strength_ratios = [], []
    try:
        if hasattr(gmm_class, 'COEFFS'):
            coeffs = gmm_class.COEFFS
            for key in coeffs.keys():
                if 'R=' in key:
                    try:
                        ratio = float(key.split('R=')[1].split(',')[0].strip())
                        if ratio not in strength_ratios: strength_ratios.append(ratio)
                    except: continue
            if strength_ratios:
                first_key = list(coeffs.keys())[0]
                coeff_table = coeffs[first_key]
                for imt_key in coeff_table.sa_coeffs.keys():
                    if hasattr(imt_key, 'period') and float(imt_key.period) not in periods: periods.append(float(imt_key.period))
        if not periods: periods = [0.04, 0.06, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
        if not strength_ratios: strength_ratios = [1.5, 2.0, 3.0, 4.0, 6.0]
    except:
        periods = [0.04, 0.06, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
        strength_ratios = [1.5, 2.0, 3.0, 4.0, 6.0]
    periods.sort(); strength_ratios.sort()
    return periods, strength_ratios

def get_sa_periods(gmm_class):
    periods = []
    try:
        for attr_name in dir(gmm_class):
            try:
                attr = getattr(gmm_class, attr_name)
                if isinstance(attr, CoeffsTable):
                    coeffs_dict = attr.sa_coeffs if hasattr(attr, 'sa_coeffs') else getattr(attr, '_coeffs', {})
                    if coeffs_dict:
                        for key in coeffs_dict.keys():
                            if isinstance(key, (int, float)) and key > 0:
                                periods.append(float(key))
                            elif hasattr(key, 'period') and key.period > 0:
                                periods.append(float(key.period))
            except: 
                continue
        if not periods:
            periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    except:
        periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    return sorted(list(set(periods)))

def test_imt_support(gmm, sctx, rctx, dctx, imt_obj):
    try:
        mean, stddevs = gmm.get_mean_and_stddevs(sctx, rctx, dctx, imt_obj, [const.StdDev.TOTAL])
        return True, mean, stddevs
    except Exception:
        return False, None, None

def create_response_spectrum_plot(sa_results, gmm_name, output_path="/tmp/response_spectrum.png"):
    if not sa_results: return None
    try:
        sa_results.sort(key=lambda r: r['sa_period'])
        periods = [r["sa_period"] for r in sa_results]
        mean_values = np.array([r["mean_value"] for r in sa_results])
        stddev_ln_values = np.array([r["stddev_total_ln"] for r in sa_results])
        plt.figure(figsize=(10, 6)); plt.fill_between(periods, mean_values * np.exp(-stddev_ln_values), mean_values * np.exp(stddev_ln_values), color='royalblue', alpha=0.2, label='Mean ± 1σ range'); plt.plot(periods, mean_values, 'b-o', color='darkblue', label='Mean SA', linewidth=2, markersize=4); plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5); plt.xlabel('Period (s)'); plt.ylabel('Spectral Acceleration (g)'); plt.title(f'Response Spectrum - {gmm_name}', fontweight='bold'); plt.legend(); plt.tight_layout(); plt.savefig(output_path, dpi=150); plt.close()
        return output_path
    except Exception:
        return None

def create_avg_sa_spectrum_plot(avg_sa_results, gmm_name, output_path="/tmp/avg_sa_spectrum.png"):
    if not avg_sa_results: return None
    try:
        avg_sa_results.sort(key=lambda r: r['sa_period'])
        periods = [r["sa_period"] for r in avg_sa_results]
        mean_values = np.array([r["mean_value"] for r in avg_sa_results])
        stddev_ln_values = np.array([r["stddev_total_ln"] for r in avg_sa_results])
        plt.figure(figsize=(10, 6))
        plt.fill_between(periods, mean_values * np.exp(-stddev_ln_values), mean_values * np.exp(stddev_ln_values), color='darkorange', alpha=0.2, label='Mean ± 1σ range')
        plt.plot(periods, mean_values, '-o', color='orangered', label='Mean AvgSA', linewidth=2, markersize=4)
        plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.xlabel('Period (s)'); plt.ylabel('Average Spectral Acceleration (g)'); plt.title(f'Average SA Spectrum - {gmm_name}', fontweight='bold'); plt.legend(); plt.tight_layout()
        plt.savefig(output_path, dpi=150); plt.close()
        return output_path
    except Exception as e:
        print(f"Erreur lors de la création du graphique AvgSA: {e}", file=sys.stderr)
        return None

def create_vh_ratio_spectrum_plot(vhr_results, gmm_name, output_path="/tmp/vh_ratio_spectrum.png"):
    if not vhr_results: return None
    try:
        vhr_results.sort(key=lambda r: r['sa_period'])
        periods = [r["sa_period"] for r in vhr_results]
        mean_values = np.array([r["mean_value"] for r in vhr_results])
        stddev_ln_values = np.array([r["stddev_total_ln"] for r in vhr_results])
    except Exception as e:
        print(f"Error preparing V/H ratio plot data: {e}", file=sys.stderr)
        return None
    
    plt.figure(figsize=(10, 6))
    plt.fill_between(periods, mean_values * np.exp(-stddev_ln_values), mean_values * np.exp(stddev_ln_values), color='mediumseagreen', alpha=0.2, label='Mean ± 1σ range')
    plt.plot(periods, mean_values, '-o', color='darkgreen', label='Mean V/H Ratio', linewidth=2, markersize=4)
    plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.xlabel('Period (s)'); plt.ylabel('V/H Ratio'); plt.title(f'V/H Ratio Spectrum - {gmm_name}', fontweight='bold')
    plt.legend(); plt.tight_layout(); plt.savefig(output_path, dpi=150); plt.close()
    return output_path

def create_sdi_spectrum_plot(sdi_results, gmm_name, output_path="/tmp/sdi_spectrum.png"):
    if not sdi_results: return None
    try:
        ratios_data = {}; colors = ['darkblue', 'darkgreen', 'darkred', 'orange', 'purple']
        for result in sdi_results:
            ratio = result["strength_ratio"]; 
            if ratio not in ratios_data: ratios_data[ratio] = []
            ratios_data[ratio].append(result)
        plt.figure(figsize=(12, 8))
        for i, ratio in enumerate(sorted(ratios_data.keys())):
            data = sorted(ratios_data[ratio], key=lambda r: r['sdi_period']); periods = [r["sdi_period"] for r in data]; mean_values = np.array([r["mean_value"] for r in data]); stddev_ln = np.array([r["stddev_total_ln"] for r in data]); color = colors[i % len(colors)]; plt.fill_between(periods, mean_values * np.exp(-stddev_ln), mean_values * np.exp(stddev_ln), color=color, alpha=0.1); plt.plot(periods, mean_values, '-o', color=color, label=f'Mean SDi (R={ratio})', linewidth=2, markersize=4)
        plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5); plt.xlabel('Period (s)'); plt.ylabel('Inelastic Spectral Displacement (cm)'); plt.title(f'Inelastic Displacement Spectra - {gmm_name}', fontweight='bold'); plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left'); plt.tight_layout(); plt.savefig(output_path, dpi=150); plt.close(); return output_path
    except Exception as e: print(f"Error creating SDi plot: {e}", file=sys.stderr); return None

def create_FAS_spectrum_plot(fas_results, gmm_name, output_path="/tmp/fas_spectrum.png"):
    if not fas_results: return None
    try:
        fas_results.sort(key=lambda r: r['frequency']); frequencies = [r["frequency"] for r in fas_results]; mean_values = np.array([r["mean_value"] for r in fas_results]); stddev_ln_values = np.array([r["stddev_total_ln"] for r in fas_results]); plt.figure(figsize=(10, 6)); plt.fill_between(frequencies, mean_values * np.exp(-stddev_ln_values), mean_values * np.exp(stddev_ln_values), color='darkgreen', alpha=0.2, label='Mean ± 1σ range'); plt.plot(frequencies, mean_values, '-o', color='darkgreen', label='Mean FAS', linewidth=2, markersize=4); plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5); plt.xlabel('Frequency (Hz)'); plt.ylabel('Fourier Amplitude Spectrum (g.s)'); plt.title(f'Fourier Amplitude Spectrum - {gmm_name}', fontweight='bold'); plt.legend(); plt.tight_layout(); plt.savefig(output_path, dpi=150); plt.close(); return output_path
    except Exception as e: print(f"Error creating FAS plot: {e}", file=sys.stderr); return None

def create_EAS_spectrum_plot(eas_results, gmm_name, output_path="/tmp/eas_spectrum.png"):
    if not eas_results: return None
    try:
        eas_results.sort(key=lambda r: r['frequency']); frequencies = [r["frequency"] for r in eas_results]; mean_values = np.array([r["mean_value"] for r in eas_results]); stddev_ln_values = np.array([r["stddev_total_ln"] for r in eas_results]); plt.figure(figsize=(10, 6)); plt.fill_between(frequencies, mean_values * np.exp(-stddev_ln_values), mean_values * np.exp(stddev_ln_values), color='purple', alpha=0.2, label='Mean ± 1σ range'); plt.plot(frequencies, mean_values, '-o', color='purple', label='Mean EAS', linewidth=2, markersize=4); plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5); plt.xlabel('Frequency (Hz)'); plt.ylabel('Effective Amplitude Spectrum (cm/s)'); plt.title(f'Effective Amplitude Spectrum - {gmm_name}', fontweight='bold'); plt.legend(); plt.tight_layout(); plt.savefig(output_path, dpi=150); plt.close(); return output_path
    except Exception as e: print(f"Error creating EAS plot: {e}", file=sys.stderr); return None

def create_DRVT_spectrum_plot(drvt_results, gmm_name, output_path="/tmp/drvt_spectrum.png"):
    if not drvt_results: return None
    try:
        drvt_results.sort(key=lambda r: r['frequency'])
        frequencies = [r["frequency"] for r in drvt_results]
        mean_values = np.array([r["mean_value"] for r in drvt_results])
        stddev_ln_values = np.array([r["stddev_total_ln"] for r in drvt_results])
        plt.figure(figsize=(10, 6))
        plt.fill_between(frequencies, np.exp(np.log(mean_values) - stddev_ln_values), np.exp(np.log(mean_values) + stddev_ln_values), color='teal', alpha=0.2, label='Mean ± 1σ range')
        plt.plot(frequencies, mean_values, '-o', color='darkcyan', label='Mean DRVT', linewidth=2, markersize=4)
        plt.xscale('log'); plt.yscale('log'); plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.xlabel('Frequency (Hz)'); plt.ylabel('Duration (s)')
        plt.title(f'Duration Spectrum (DRVT) - {gmm_name}', fontweight='bold'); plt.legend(); plt.tight_layout()
        plt.savefig(output_path, dpi=150); plt.close()
        return output_path
    except Exception as e:
        print(f"Erreur lors de la préparation des données pour le tracé DRVT: {e}", file=sys.stderr)
        return None


def run_gmm_calculation(args):
    try:
        keys = [k.strip() for k in args.keys.split(',')]
        raw_values = [v.strip() for v in args.values.split(',')]
        values = []
        for v in raw_values:
            try: values.append(float(v))
            except ValueError: values.append(v)
        user_params = dict(zip(keys, values))
        full_module_name = f"openquake.hazardlib.gsim.{args.module}"
        module = importlib.import_module(full_module_name)
        gmm_class = getattr(module, args.class_name)

        # =================================================================================
        # ### DÉBUT DE LA MODIFICATION V8.5 : MAPPING NUMÉRIQUE POUR LA GÉOLOGIE ###
        # =================================================================================
        # On applique cette logique UNIQUEMENT pour la classe GMM spécifiée.
        if args.class_name == 'Weatherill2024ESHM20SlopeGeologyAvgSA':
            # Dictionnaire de correspondance : l'utilisateur entre le chiffre, le script utilise le texte.
            # L'utilisateur entrera un nombre de 1 à 8.
            geology_map = (
                "CENOZOIC", "HOLOCENE", "JURASSIC-TRIASSIC", "CRETACEOUS",
                "PALEOZOIC", "PLEISTOCENE", "PRECAMBRIAN", "UNKNOWN"
            )
            
            # On vérifie si le paramètre 'geology' a été fourni.
            if 'geology' in user_params:
                try:
                    # On convertit l'entrée de l'utilisateur en index (ex: 1 -> index 0)
                    geology_index = int(user_params['geology']) - 1
                    
                    # On vérifie si l'index est valide.
                    if 0 <= geology_index < len(geology_map):
                        # On remplace le nombre par la chaîne de caractères correspondante.
                        user_params['geology'] = geology_map[geology_index]
                    else:
                        # Si le nombre est hors limites, on lève une erreur.
                        raise ValueError(f"Numéro de géologie '{user_params['geology']}' invalide. Doit être entre 1 et {len(geology_map)}.")
                except (ValueError, TypeError):
                    # Si l'entrée n'est pas un nombre, on lève une erreur.
                    raise ValueError(f"L'entrée pour la géologie '{user_params['geology']}' est invalide. Un nombre entier est attendu.")
        # =================================================================================
        # ### FIN DE LA MODIFICATION V8.5 ###
        # =================================================================================


        constructor_args = {}
        # Gestion spécifique pour HassaniAtkinson2018
        if 'HassaniAtkinson2018' in args.class_name:
            # Paramètres obligatoires avec valeurs par défaut
            constructor_args['d_sigma'] = user_params.pop('d_sigma')  # stress drop en bars
            constructor_args['kappa0'] = user_params.pop('kappa0',)    # kappa0 en secondes
          
            # Paramètre optionnel pour l'atténuation anélastique
            if 'gamma_fle' in user_params:
                constructor_args['gamma_fle'] = user_params.pop('gamma_fle')
        elif 'MacedoEtAl2019' in args.class_name:
            if 'region' in user_params: constructor_args['region'] = str(user_params.pop('region'))
            constructor_args['gmpe'] = {'AbrahamsonEtAl2015SInter': {}}

        # Gestion des paramètres spécifiques pour différents GMMs
        # Logique pour les paramètres region et siteclass
        if ('KothaEtAl2020ESHM20' in args.class_name or 'Weatherill2024ESHM20' in args.class_name):
            if 'region' not in user_params:
                user_params['region'] = 3
            else:
                # S'assurer que region est un entier pour ces modèles
                try:
                    user_params['region'] = int(float(user_params['region']))
                except (ValueError, TypeError):
                    user_params['region'] = 3
        
        if 'LanzanoEtAl2020_Cluster' in args.class_name:
            if 'siteclass' not in user_params:
                user_params['siteclass'] = "3"
            # S'assurer que siteclass est une string pour ce modèle
            user_params['siteclass'] = str(user_params['siteclass'])


        init_params = inspect.signature(gmm_class.__init__).parameters
        if 'gmpe_name' in init_params:
            constructor_args['gmpe_name'] = 'AbrahamsonSilva2008'
        gmm = gmm_class(**constructor_args)
    except Exception as e:
        print(json.dumps({"success": False, "error": f"GMM Import Error: {e}"}), file=sys.stderr)
        sys.exit(1)
    
    is_vh_ratio_model = (hasattr(gmm_class, 'DEFINED_FOR_INTENSITY_MEASURE_COMPONENT') and
                       gmm_class.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT == const.IMC.VERTICAL_TO_HORIZONTAL_RATIO)
    try:
        if ('KothaEtAl2020ESHM20' in args.class_name or 'Weatherill2024ESHM20' in args.class_name) and 'region' in user_params:
            try: user_params['region'] = int(float(user_params['region']))
            except (ValueError, TypeError): user_params['region'] = 3

        req_sites = getattr(gmm_class, 'REQUIRES_SITES_PARAMETERS', set())
        req_rupture = getattr(gmm_class, 'REQUIRES_RUPTURE_PARAMETERS', set())
        req_distances = getattr(gmm_class, 'REQUIRES_DISTANCES', set())
        all_required = req_sites | req_rupture | req_distances
        missing = [p for p in all_required if p not in user_params]
        if missing:
            print(json.dumps({"success": False, "error": f"Missing parameters: {', '.join(missing)}"}), file=sys.stderr)
            sys.exit(1)
        
        sctx, rctx, dctx = SitesContext(), RuptureContext(), DistancesContext()
        for p in req_sites: setattr(sctx, p, np.array([user_params[p]]))
        for p in req_rupture: setattr(rctx, p, np.array([user_params[p]]))
        for p in req_distances: setattr(dctx, p, np.array([user_params[p]]))
        sctx.sids = np.array([1])
    except Exception as e:
        print(json.dumps({"success": False, "error": f"Context Error: {e}"}), file=sys.stderr)
        sys.exit(1)

    try:
        supported_imts = get_supported_imts(gmm_class)
    except Exception as e:
        print(json.dumps({"success": False, "error": f"IMT Detection Error: {e}"}), file=sys.stderr)
        sys.exit(1)
        
    results, successful_imts, failed_imts = [], [], []
    for imt_name, imt_params in supported_imts:
        display_name = ""
        try:
            if imt_name == "SA": imt_obj, display_name = imt.SA(imt_params), f"SA({imt_params}s)"
            elif imt_name == "AvgSA": imt_obj, display_name = imt.AvgSA(imt_params), f"AvgSA({imt_params}s)"
            elif imt_name == "FAS": imt_obj, display_name = imt.FAS(imt_params), f"FAS({imt_params}Hz)"
            elif imt_name == "EAS": imt_obj, display_name = imt.EAS(imt_params), f"EAS({imt_params}Hz)"
            elif imt_name == "DRVT": imt_obj, display_name = imt.DRVT(imt_params), f"DRVT({imt_params}Hz)"
            elif imt_name == "SDi": imt_obj, display_name = imt.SDi(imt_params["period"], imt_params["strength_ratio"]), f"SDi({imt_params['period']}s, R={imt_params['strength_ratio']})"
            else: imt_obj, display_name = getattr(imt, imt_name)(), imt_name
            
            is_supported, mean, stddevs = test_imt_support(gmm, sctx, rctx, dctx, imt_obj)
            
            stddev_val = float(stddevs[0][0]) if stddevs and len(stddevs) > 0 and len(stddevs[0]) > 0 else 9999.0
            
            if is_supported and stddev_val < 9000:
                if is_vh_ratio_model:
                    mean_value = float(np.exp(mean[0]))
                    unit = ""
                    imt_name_orig = imt_name
                    imt_name = f"VHR_{imt_name_orig}"
                    display_name = f"VHR({imt_params}s)" if imt_name_orig in ["SA", "AvgSA"] else f"VHR_{display_name}"
                else:
                    # ### DÉBUT DU BLOC CORRIGÉ (INDENTATION PROPRE) ###
                    if "MorikawaFujiwara2013" in args.class_name and imt_name in ["JMA", "PGV"]:
                        mean_value = float(np.log10(np.exp(mean[0]) * 980.665))
                    elif imt_name in ["MMI", "JMA"]:
                        mean_value = float(mean[0])
                    elif 'ConvertitoEtAl2012' in args.class_name or 'TusaLangerAzzaro2019' in args.class_name:
                        mean_value = float(np.exp(mean[0]) / 1000.0)
                    else:
                        mean_value = float(np.exp(mean[0]))
                    # ### FIN DU BLOC CORRIGÉ ###
                    unit_map = {"PGA":"g", "PGV":"cm/s", "PGD":"cm", "IA":"m/s", "CAV":"g⋅s", "SA":"g", "AvgSA": "g", "FAS":"g.s", "SDi":"cm", "EAS":"cm/s", "DRVT":"s"}
                    duration_units = {"RSD575", "RSD595", "RSD2080", "D5_75", "D5_95"}
                    unit = unit_map.get(imt_name, "s" if imt_name in duration_units or imt_name == "DRVT" else "")
                
                result_data = {"imt": imt_name, "display_name": display_name, "mean_ln": round(float(mean[0]), 4),
                               "mean_value": round(mean_value, 6), "unit": unit, "stddev_total_ln": round(stddev_val, 4), "success": True}
                
                if imt_name in ["SA", "VHR_SA", "AvgSA"]: result_data["sa_period"] = imt_params
                elif imt_name in ["FAS", "EAS", "DRVT"]: result_data["frequency"] = imt_params
                elif imt_name == "SDi": result_data.update({"sdi_period": imt_params["period"], "strength_ratio": imt_params["strength_ratio"]})
                results.append(result_data)
                successful_imts.append(display_name)
            else:
                failed_imts.append(display_name)
        except Exception as e:
            if not display_name: display_name = f"{imt_name}({imt_params})"
            failed_imts.append(f"{display_name} (Error: {type(e).__name__})"); continue
            
    def sort_key(r):
        order = {"PGA":0, "PGV":1, "PGD":2, "LSD":3, "MMI": 3.5, "JMA": 3.6, "VHR_PGA": 3.7, "VHR_PGV": 3.8}
        if r["imt"] in order: return (order[r["imt"]], 0)
        if r["imt"] in ["SA", "AvgSA"]: return (6, r.get("sa_period", 0))
        if r["imt"] in ["FAS", "EAS", "DRVT"]: return (7, r.get("frequency", 0))
        return (9, r["imt"])
    results.sort(key=sort_key)
    
    sa_results = [r for r in results if r["imt"] == "SA"]
    avg_sa_results = [r for r in results if r["imt"] == "AvgSA"]
    vhr_sa_results = [r for r in results if r["imt"] == "VHR_SA"] 
    sdi_results = [r for r in results if r["imt"] == "SDi"]
    fas_results = [r for r in results if r["imt"] == "FAS"]
    eas_results = [r for r in results if r["imt"] == "EAS"]
    drvt_results = [r for r in results if r["imt"] == "DRVT"]
    plot_path = None
    
    if sa_results and len(sa_results) >= 3:
        plot_path = create_response_spectrum_plot(sa_results, args.class_name)
    elif avg_sa_results and len(avg_sa_results) >= 3:
        plot_path = create_avg_sa_spectrum_plot(avg_sa_results, args.class_name)
    elif vhr_sa_results and len(vhr_sa_results) >= 3:
        plot_path = create_vh_ratio_spectrum_plot(vhr_sa_results, args.class_name)
    elif drvt_results and len(drvt_results) >= 3:
        plot_path = create_DRVT_spectrum_plot(drvt_results, args.class_name)
    elif fas_results and len(fas_results) >= 3:
        plot_path = create_FAS_spectrum_plot(fas_results, args.class_name)
    elif eas_results and len(eas_results) >= 3:
        plot_path = create_EAS_spectrum_plot(eas_results, args.class_name)
    elif sdi_results and len(sdi_results) >= 3:
        plot_path = create_sdi_spectrum_plot(sdi_results, args.class_name)
    
    response = {
        "success": True, "gmm": args.class_name, "module": args.module,
        "total_imts_tested": len(supported_imts), "successful_imts_count": len(results),
        "failed_imts_count": len(failed_imts), "successful_imts": successful_imts,
        "failed_imts": failed_imts or None, "imt_results": results,
        "sa_count": len(sa_results), "avg_sa_count": len(avg_sa_results),
        "vhr_sa_count": len(vhr_sa_results), "sdi_count": len(sdi_results),
        "fas_count": len(fas_results), "eas_count": len(eas_results),
        "drvt_count": len(drvt_results), "plot_path": plot_path
    }
    print(json.dumps(response, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calcul dynamique et intelligent de tous les IMTs supportés par un GMM OpenQuake")
    parser.add_argument('--module', required=True, help="Module GMM")
    parser.add_argument('--class_name', required=True, help="Nom de la classe GMM")
    parser.add_argument('--keys', required=True, help="Clés des paramètres séparées par virgules")
    parser.add_argument('--values', required=True, help="Valeurs des paramètres séparées par virgules")
    args = parser.parse_args()
    run_gmm_calculation(args)