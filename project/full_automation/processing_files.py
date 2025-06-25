import os
import glob
import time
from collections import defaultdict

from parsing_stream import parsing_stream
from parsing_err_file import parse_err
from parsing_UC_files import parse_UC_file
from resolution_cutoff_determination import calculating_max_res_from_Rsplit_CCstar_dat
from preparation_for_statistics_calculations import get_UC
from get_SNR_resolution import get_d_at_snr_one
from get_CC_threshold_resolution import get_d_at_cc_threshold

indexes = [
    'Num. patterns/hits', 'Indexed patterns/crystals', 'Resolution', 'Rsplit(%)',
    'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness(%)', 'Multiplicity',
    'Total Measurements', 'Unique Reflections', 'Wilson B-factor',
    'Resolution SNR=1', 'Resolution CC>=0.3', 'a,b,c,alpha,betta,gamma'
]

def processing_statistics_for_run(
    name_of_run, data_info_for_the_current_run, hkl_file, main_path,
    is_extended=False, cell_path=None, is_refining=False
):
    print(f'Processing {name_of_run}')
    data_info = defaultdict(dict)
    data_info[name_of_run] = {key: '' for key in indexes}

    base, _ = os.path.splitext(hkl_file)
    CCstar_dat_file = f"{base}_CCstar.dat"
    Rsplit_dat_file = f"{base}_Rsplit.dat"
    SNR_dat_file = f"{base}_SNR.dat"
    CC_dat_file = f"{base}_CC.dat"
    stream_file = f"{base}.stream" if "_offset_" not in base else base.split("_offset_")[0] + '.stream'

    # Locate UC file
    UC_file = get_UC(hkl_file)
    if UC_file and not os.path.exists(UC_file):
        alt_path = os.path.join(os.path.dirname(hkl_file), os.path.basename(UC_file))
        UC_file = alt_path if os.path.exists(alt_path) else None

    if cell_path:
        base_name = os.path.basename(hkl_file)
        for ext in ['cell', 'pdb']:
            candidate = os.path.join(cell_path, base_name.replace("hkl", ext))
            if os.path.exists(candidate):
                UC_file = candidate
                break

    # Parse stream file
    chunks, hits, indexed_patterns, indexed = parsing_stream(stream_file)
    data_info[name_of_run]['Num. patterns/hits'] = f"{chunks}/{hits}"
    data_info[name_of_run]['Indexed patterns/crystals'] = f"{indexed_patterns}/{indexed}"

    # Parse unit cell
    a, b, c, al, be, ga = parse_UC_file(UC_file) if UC_file else (None,) * 6
    data_info[name_of_run]['a,b,c,alpha,betta,gamma'] = (a, b, c, al, be, ga)

    # Wait for required files
    for file_path in [Rsplit_dat_file, SNR_dat_file, CC_dat_file]:
        while not os.path.exists(file_path) or os.stat(file_path).st_size == 0:
            time.sleep(5)

    # Resolution metrics
    d_snr = get_d_at_snr_one(SNR_dat_file)
    d_cc = get_d_at_cc_threshold(CC_dat_file)
    max_res = calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)

    Rwork = data_info_for_the_current_run[name_of_run].get('Rwork', '')
    Rfree = data_info_for_the_current_run[name_of_run].get('Rfree', '')
    resolution_cut_off_high = data_info_for_the_current_run[name_of_run].get('resolution_cut_off_high', '')
    resolution_cut_off_low = data_info_for_the_current_run[name_of_run].get('resolution_low', '')

    if Rwork is None and is_refining:
        print('Calling refinement...')
        # Call refinement process here
        pass

    data_info[name_of_run].update({
        'Rwork/Rfree': f'{Rwork}/{Rfree}',
        'Refinement resolution cut-off high': str(resolution_cut_off_high),
        'Refinement resolution cut-off low': str(resolution_cut_off_low),
        'Resolution SNR=1': str(d_snr),
        'Resolution CC>=0.3': str(d_cc),
        'CC* intersects with Rsplit at': str(max_res)
    })

    # Parse error file
    err_pattern = os.path.join(
        main_path, '**', os.path.basename(CCstar_dat_file).replace("_CCstar.dat", "*.err")
    )
    err_files = glob.glob(err_pattern, recursive=True)
    if err_files:
        latest_err = max(err_files, key=os.path.getctime)
        data_info = parse_err(data_info, name_of_run, latest_err, CCstar_dat_file)

    # Convert to MTZ if UC file exists
    mtz_file = hkl_file.replace("hkl", "mtz")
    if UC_file:
        mtz_command = f'get_hkl -i {hkl_file} -p {UC_file} -o {mtz_file} --output-format=mtz | tee conversion_log.txt'
        os.system(mtz_command)
    else:
        mtz_command = ''

    # Extended data
    if is_extended:
        data_info[name_of_run].update({
            'UC_file': UC_file,
            'N_patterns': str(chunks),
            'N_hits': str(hits),
            'Indexed_patterns': str(indexed_patterns),
            'Indexed_crystals': str(indexed),
            'a': a, 'b': b, 'c': c,
            'alpha': al, 'betta': be, 'gamma': ga,
            'mtz_command': mtz_command,
            'mtz_with_the_defined_resolution': f"{mtz_file} {d_snr}" if d_snr else ''
        })

    return data_info
