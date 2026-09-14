#!/usr/bin/env python

'''A script to combine the toys run by TRExFitter and compute the limits.

Typical usage example:

    python combine_toys -i /path/to/inputFile.txt -o /path/to/outputFile.root -n result_mu

'''

__author__ = "Jackson Carl Burzynski"
__email__ = "jackson.carl.burzynski@cern.ch"

import ROOT
import argparse
import array

def combine_toys(files, result_name, write_output=False, output_file="combined_results.root"):
    '''
    Returns limits computed from combining toys from the provided input file.

    Retrieves the HypoTestResult object from each input file and combines them
    using HypoTestInverterResult::Add. Returns the expected/observed limits as
    a dictionary and optionally writes a TTree with the values.

            Args:
                    files: A List of TFiles to process
                    result_name: The name of the HypoTestResult object to use
                    write_output: Option to write out a TTree
                    output_file: Name of the output file to write

            Returns:
                    A dict with the computed limits. For example:

                    {'p2':  0.0194, 
                     'p1':  0.0169, 
                     'exp': 0.0153, 
                     'm1':  0.0129, 
                     'm2':  0.0115, 
                     'obs': 0.0156}
    '''

    limits = {'exp': None, 'p1': None, 'p2': None, 'm1': None, 'm2': None, 'obs': None}
    hypo_test_result = None

    # add results together using HypoTestInverterResult::Add
    for f in files:
        if not f:
            continue
        if not hypo_test_result:
            hypo_test_result = f.Get(result_name)
        else:
            hypo_test_result.Add(f.Get(result_name))

    limits['obs'] = hypo_test_result.UpperLimit()
    limits['m2']  = hypo_test_result.GetExpectedUpperLimit(-2)
    limits['m1']  = hypo_test_result.GetExpectedUpperLimit(-1)
    limits['exp'] = hypo_test_result.GetExpectedUpperLimit(0)
    limits['p1']  = hypo_test_result.GetExpectedUpperLimit(1)
    limits['p2']  = hypo_test_result.GetExpectedUpperLimit(2)
       
    if write_output:
        outfile = ROOT.TFile.Open(output_file, "RECREATE")

        c1 = ROOT.TCanvas("c1")
        c1.SetLogy(False)

        plot = ROOT.RooStats.HypoTestInverterPlot("HTI_Result_Plot","test",hypo_test_result);
        plot.Draw("")
        c1.SaveAs("brasilianFlag.pdf")

        tree = ROOT.TTree("stats","")

        obs  = array.array('d', [limits['obs']])
        m2   = array.array('d', [limits['m2']])
        m1   = array.array('d', [limits['m1']])
        exp  = array.array('d', [limits['exp']])
        p1   = array.array('d', [limits['p1']])
        p2   = array.array('d', [limits['p2']])

        tree.Branch("exp_upperlimit", exp, "exp_upperlimit/D")
        tree.Branch("obs_upperlimit", obs, "obs_upperlimit/D")
        tree.Branch("exp_upperlimit_plus1", p1, "exp_upperlimit_plus1/D") 
        tree.Branch("exp_upperlimit_plus2", p2, "exp_upperlimit_plus2/D")
        tree.Branch("exp_upperlimit_minus1", m1, "exp_upperlimit_minus1/D")
        tree.Branch("exp_upperlimit_minus2", m2, "exp_upperlimit_minus2/D")

        tree.Fill()
        tree.Write()
        hypo_test_result.Write()
        outfile.Close()
    return limits

if __name__ == "__main__":

    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Input file (.txt) with list of files to combine')
    parser.add_argument('-o', '--output_file', help='Output file (.root) to \
                                           store the limit results computed from the combined toys')
    parser.add_argument('-n', '--result_name', help='Name of the HypoTestResult object (e.g. "result_mu")')
    args = parser.parse_args()

    file_list = []

    # create a list of the input file names to pass to combine_toys
    with open(args.input_file) as tmp:
        for file_name in tmp:
            try:
                f = ROOT.TFile(file_name.strip(), "READ")
            except OSError as error:
                print('Error in opening ROOT file: {}'.format(file_name.strip()))
                continue

            file_list.append(f)

    result = combine_toys(file_list, args.result_name, write_output=True, output_file=args.output_file)
    print("Combined limit results:")
    print("\t Expected: {}".format(result['exp']))
    print("\t Expected + 2sigma: {}".format(result['p2']))
    print("\t Expected + 1sigma: {}".format(result['p1']))
    print("\t Expected - 1sigma: {}".format(result['m1']))
    print("\t Expected - 2sigma: {}".format(result['m2']))
    print("\t Observed: {}".format(result['obs']))


