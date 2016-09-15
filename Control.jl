# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: Control
# Description: Calculate the transcriptional control array at time t
# Generated on: 2016-09-14T20:16:02
#
# Input arguments:
# t::Float64 => Current time value (scalar) 
# x::Array{Float64,1} => State array (number_of_species x 1) 
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model parameters 
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t 
# ----------------------------------------------------------------------------------- #
function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})

	# initialize the control - 
	control_array = zeros(17)

	# Alias the species - 
	gene_AP1 = x[1]
	gene_AhR = x[2]
	gene_CD11b = x[3]
	gene_CD14 = x[4]
	gene_CD38 = x[5]
	gene_CEBPa = x[6]
	gene_E2F = x[7]
	gene_EGR1 = x[8]
	gene_GFI1 = x[9]
	gene_IRF1 = x[10]
	gene_OCT1 = x[11]
	gene_OCT4 = x[12]
	gene_P21 = x[13]
	gene_P47Phox = x[14]
	gene_PPARg = x[15]
	gene_PU1 = x[16]
	gene_RARa = x[17]
	mRNA_gene_AP1 = x[18]
	mRNA_gene_AhR = x[19]
	mRNA_gene_CD11b = x[20]
	mRNA_gene_CD14 = x[21]
	mRNA_gene_CD38 = x[22]
	mRNA_gene_CEBPa = x[23]
	mRNA_gene_E2F = x[24]
	mRNA_gene_EGR1 = x[25]
	mRNA_gene_GFI1 = x[26]
	mRNA_gene_IRF1 = x[27]
	mRNA_gene_OCT1 = x[28]
	mRNA_gene_OCT4 = x[29]
	mRNA_gene_P21 = x[30]
	mRNA_gene_P47Phox = x[31]
	mRNA_gene_PPARg = x[32]
	mRNA_gene_PU1 = x[33]
	mRNA_gene_RARa = x[34]
	protein_gene_AP1 = x[35]
	protein_gene_AhR = x[36]
	protein_gene_CD11b = x[37]
	protein_gene_CD14 = x[38]
	protein_gene_CD38 = x[39]
	protein_gene_CEBPa = x[40]
	protein_gene_E2F = x[41]
	protein_gene_EGR1 = x[42]
	protein_gene_GFI1 = x[43]
	protein_gene_IRF1 = x[44]
	protein_gene_OCT1 = x[45]
	protein_gene_OCT4 = x[46]
	protein_gene_P21 = x[47]
	protein_gene_P47Phox = x[48]
	protein_gene_PPARg = x[49]
	protein_gene_PU1 = x[50]
	protein_gene_RARa = x[51]

	# Alias the binding parameters - 
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_gene_AP1_gene_AhR = binding_parameter_dictionary["n_gene_AP1_gene_AhR"]
	K_gene_AP1_gene_AhR = binding_parameter_dictionary["K_gene_AP1_gene_AhR"]
	n_gene_AP1_gene_PU1 = binding_parameter_dictionary["n_gene_AP1_gene_PU1"]
	K_gene_AP1_gene_PU1 = binding_parameter_dictionary["K_gene_AP1_gene_PU1"]
	n_gene_AP1_gene_PPARg = binding_parameter_dictionary["n_gene_AP1_gene_PPARg"]
	K_gene_AP1_gene_PPARg = binding_parameter_dictionary["K_gene_AP1_gene_PPARg"]
	n_gene_AhR_gene_RARa = binding_parameter_dictionary["n_gene_AhR_gene_RARa"]
	K_gene_AhR_gene_RARa = binding_parameter_dictionary["K_gene_AhR_gene_RARa"]
	n_gene_CD11b_gene_PU1 = binding_parameter_dictionary["n_gene_CD11b_gene_PU1"]
	K_gene_CD11b_gene_PU1 = binding_parameter_dictionary["K_gene_CD11b_gene_PU1"]
	n_gene_CD14_gene_EGR1 = binding_parameter_dictionary["n_gene_CD14_gene_EGR1"]
	K_gene_CD14_gene_EGR1 = binding_parameter_dictionary["K_gene_CD14_gene_EGR1"]
	n_gene_CD14_gene_CEBPa = binding_parameter_dictionary["n_gene_CD14_gene_CEBPa"]
	K_gene_CD14_gene_CEBPa = binding_parameter_dictionary["K_gene_CD14_gene_CEBPa"]
	n_gene_CD14_gene_PPARg = binding_parameter_dictionary["n_gene_CD14_gene_PPARg"]
	K_gene_CD14_gene_PPARg = binding_parameter_dictionary["K_gene_CD14_gene_PPARg"]
	n_gene_CD38_gene_RARa = binding_parameter_dictionary["n_gene_CD38_gene_RARa"]
	K_gene_CD38_gene_RARa = binding_parameter_dictionary["K_gene_CD38_gene_RARa"]
	n_gene_CD38_gene_PPARg = binding_parameter_dictionary["n_gene_CD38_gene_PPARg"]
	K_gene_CD38_gene_PPARg = binding_parameter_dictionary["K_gene_CD38_gene_PPARg"]
	n_gene_CD38_gene_IRF1 = binding_parameter_dictionary["n_gene_CD38_gene_IRF1"]
	K_gene_CD38_gene_IRF1 = binding_parameter_dictionary["K_gene_CD38_gene_IRF1"]
	n_gene_CEBPa_gene_RARa = binding_parameter_dictionary["n_gene_CEBPa_gene_RARa"]
	K_gene_CEBPa_gene_RARa = binding_parameter_dictionary["K_gene_CEBPa_gene_RARa"]
	n_gene_CEBPa_gene_PPARg = binding_parameter_dictionary["n_gene_CEBPa_gene_PPARg"]
	K_gene_CEBPa_gene_PPARg = binding_parameter_dictionary["K_gene_CEBPa_gene_PPARg"]
	n_gene_CEBPa_gene_CEBPa = binding_parameter_dictionary["n_gene_CEBPa_gene_CEBPa"]
	K_gene_CEBPa_gene_CEBPa = binding_parameter_dictionary["K_gene_CEBPa_gene_CEBPa"]
	n_gene_CEBPa_gene_GFI1 = binding_parameter_dictionary["n_gene_CEBPa_gene_GFI1"]
	K_gene_CEBPa_gene_GFI1 = binding_parameter_dictionary["K_gene_CEBPa_gene_GFI1"]
	n_gene_E2F_gene_E2F = binding_parameter_dictionary["n_gene_E2F_gene_E2F"]
	K_gene_E2F_gene_E2F = binding_parameter_dictionary["K_gene_E2F_gene_E2F"]
	n_gene_E2F_gene_PPARg = binding_parameter_dictionary["n_gene_E2F_gene_PPARg"]
	K_gene_E2F_gene_PPARg = binding_parameter_dictionary["K_gene_E2F_gene_PPARg"]
	n_gene_E2F_gene_CEBPa = binding_parameter_dictionary["n_gene_E2F_gene_CEBPa"]
	K_gene_E2F_gene_CEBPa = binding_parameter_dictionary["K_gene_E2F_gene_CEBPa"]
	n_gene_E2F_gene_GFI1 = binding_parameter_dictionary["n_gene_E2F_gene_GFI1"]
	K_gene_E2F_gene_GFI1 = binding_parameter_dictionary["K_gene_E2F_gene_GFI1"]
	n_gene_EGR1_gene_RARa = binding_parameter_dictionary["n_gene_EGR1_gene_RARa"]
	K_gene_EGR1_gene_RARa = binding_parameter_dictionary["K_gene_EGR1_gene_RARa"]
	n_gene_EGR1_gene_PU1 = binding_parameter_dictionary["n_gene_EGR1_gene_PU1"]
	K_gene_EGR1_gene_PU1 = binding_parameter_dictionary["K_gene_EGR1_gene_PU1"]
	n_gene_EGR1_gene_PPARg = binding_parameter_dictionary["n_gene_EGR1_gene_PPARg"]
	K_gene_EGR1_gene_PPARg = binding_parameter_dictionary["K_gene_EGR1_gene_PPARg"]
	n_gene_EGR1_gene_GFI1 = binding_parameter_dictionary["n_gene_EGR1_gene_GFI1"]
	K_gene_EGR1_gene_GFI1 = binding_parameter_dictionary["K_gene_EGR1_gene_GFI1"]
	n_gene_GFI1_gene_CEBPa = binding_parameter_dictionary["n_gene_GFI1_gene_CEBPa"]
	K_gene_GFI1_gene_CEBPa = binding_parameter_dictionary["K_gene_GFI1_gene_CEBPa"]
	n_gene_GFI1_gene_EGR1 = binding_parameter_dictionary["n_gene_GFI1_gene_EGR1"]
	K_gene_GFI1_gene_EGR1 = binding_parameter_dictionary["K_gene_GFI1_gene_EGR1"]
	n_gene_IRF1_gene_RARa = binding_parameter_dictionary["n_gene_IRF1_gene_RARa"]
	K_gene_IRF1_gene_RARa = binding_parameter_dictionary["K_gene_IRF1_gene_RARa"]
	n_gene_IRF1_gene_AhR = binding_parameter_dictionary["n_gene_IRF1_gene_AhR"]
	K_gene_IRF1_gene_AhR = binding_parameter_dictionary["K_gene_IRF1_gene_AhR"]
	n_gene_IRF1_gene_PPARg = binding_parameter_dictionary["n_gene_IRF1_gene_PPARg"]
	K_gene_IRF1_gene_PPARg = binding_parameter_dictionary["K_gene_IRF1_gene_PPARg"]
	n_gene_OCT1_gene_PPARg = binding_parameter_dictionary["n_gene_OCT1_gene_PPARg"]
	K_gene_OCT1_gene_PPARg = binding_parameter_dictionary["K_gene_OCT1_gene_PPARg"]
	n_gene_OCT4_gene_RARa = binding_parameter_dictionary["n_gene_OCT4_gene_RARa"]
	K_gene_OCT4_gene_RARa = binding_parameter_dictionary["K_gene_OCT4_gene_RARa"]
	n_gene_OCT4_gene_AhR = binding_parameter_dictionary["n_gene_OCT4_gene_AhR"]
	K_gene_OCT4_gene_AhR = binding_parameter_dictionary["K_gene_OCT4_gene_AhR"]
	n_gene_P21_gene_RARa = binding_parameter_dictionary["n_gene_P21_gene_RARa"]
	K_gene_P21_gene_RARa = binding_parameter_dictionary["K_gene_P21_gene_RARa"]
	n_gene_P21_gene_PPARg = binding_parameter_dictionary["n_gene_P21_gene_PPARg"]
	K_gene_P21_gene_PPARg = binding_parameter_dictionary["K_gene_P21_gene_PPARg"]
	n_gene_P21_gene_PU1 = binding_parameter_dictionary["n_gene_P21_gene_PU1"]
	K_gene_P21_gene_PU1 = binding_parameter_dictionary["K_gene_P21_gene_PU1"]
	n_gene_P21_gene_IRF1 = binding_parameter_dictionary["n_gene_P21_gene_IRF1"]
	K_gene_P21_gene_IRF1 = binding_parameter_dictionary["K_gene_P21_gene_IRF1"]
	n_gene_P21_gene_CEBPa = binding_parameter_dictionary["n_gene_P21_gene_CEBPa"]
	K_gene_P21_gene_CEBPa = binding_parameter_dictionary["K_gene_P21_gene_CEBPa"]
	n_gene_P21_gene_AP1 = binding_parameter_dictionary["n_gene_P21_gene_AP1"]
	K_gene_P21_gene_AP1 = binding_parameter_dictionary["K_gene_P21_gene_AP1"]
	n_gene_P21_gene_GFI1 = binding_parameter_dictionary["n_gene_P21_gene_GFI1"]
	K_gene_P21_gene_GFI1 = binding_parameter_dictionary["K_gene_P21_gene_GFI1"]
	n_gene_P47Phox_gene_PU1 = binding_parameter_dictionary["n_gene_P47Phox_gene_PU1"]
	K_gene_P47Phox_gene_PU1 = binding_parameter_dictionary["K_gene_P47Phox_gene_PU1"]
	n_gene_P47Phox_gene_CEBPa = binding_parameter_dictionary["n_gene_P47Phox_gene_CEBPa"]
	K_gene_P47Phox_gene_CEBPa = binding_parameter_dictionary["K_gene_P47Phox_gene_CEBPa"]
	n_gene_P47Phox_gene_PPARg = binding_parameter_dictionary["n_gene_P47Phox_gene_PPARg"]
	K_gene_P47Phox_gene_PPARg = binding_parameter_dictionary["K_gene_P47Phox_gene_PPARg"]
	n_gene_PPARg_gene_CEBPa = binding_parameter_dictionary["n_gene_PPARg_gene_CEBPa"]
	K_gene_PPARg_gene_CEBPa = binding_parameter_dictionary["K_gene_PPARg_gene_CEBPa"]
	n_gene_PPARg_gene_EGR1 = binding_parameter_dictionary["n_gene_PPARg_gene_EGR1"]
	K_gene_PPARg_gene_EGR1 = binding_parameter_dictionary["K_gene_PPARg_gene_EGR1"]
	n_gene_PPARg_gene_PU1 = binding_parameter_dictionary["n_gene_PPARg_gene_PU1"]
	K_gene_PPARg_gene_PU1 = binding_parameter_dictionary["K_gene_PPARg_gene_PU1"]
	n_gene_PPARg_gene_AP1 = binding_parameter_dictionary["n_gene_PPARg_gene_AP1"]
	K_gene_PPARg_gene_AP1 = binding_parameter_dictionary["K_gene_PPARg_gene_AP1"]
	n_gene_PU1_gene_RARa = binding_parameter_dictionary["n_gene_PU1_gene_RARa"]
	K_gene_PU1_gene_RARa = binding_parameter_dictionary["K_gene_PU1_gene_RARa"]
	n_gene_PU1_gene_AP1 = binding_parameter_dictionary["n_gene_PU1_gene_AP1"]
	K_gene_PU1_gene_AP1 = binding_parameter_dictionary["K_gene_PU1_gene_AP1"]
	n_gene_PU1_gene_OCT1 = binding_parameter_dictionary["n_gene_PU1_gene_OCT1"]
	K_gene_PU1_gene_OCT1 = binding_parameter_dictionary["K_gene_PU1_gene_OCT1"]
	n_gene_PU1_gene_AhR = binding_parameter_dictionary["n_gene_PU1_gene_AhR"]
	K_gene_PU1_gene_AhR = binding_parameter_dictionary["K_gene_PU1_gene_AhR"]
	n_gene_PU1_gene_GFI1 = binding_parameter_dictionary["n_gene_PU1_gene_GFI1"]
	K_gene_PU1_gene_GFI1 = binding_parameter_dictionary["K_gene_PU1_gene_GFI1"]
	n_gene_RARa_gene_RARa = binding_parameter_dictionary["n_gene_RARa_gene_RARa"]
	K_gene_RARa_gene_RARa = binding_parameter_dictionary["K_gene_RARa_gene_RARa"]

	# Alias the control function parameters - 
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_gene_AP1_RNAP = control_parameter_dictionary["W_gene_AP1_RNAP"]
	W_gene_AP1_gene_AhR = control_parameter_dictionary["W_gene_AP1_gene_AhR"]
	W_gene_AP1_gene_PU1 = control_parameter_dictionary["W_gene_AP1_gene_PU1"]
	W_gene_AP1_gene_PPARg = control_parameter_dictionary["W_gene_AP1_gene_PPARg"]
	W_gene_AhR_RNAP = control_parameter_dictionary["W_gene_AhR_RNAP"]
	W_gene_AhR_gene_RARa = control_parameter_dictionary["W_gene_AhR_gene_RARa"]
	W_gene_CD11b_RNAP = control_parameter_dictionary["W_gene_CD11b_RNAP"]
	W_gene_CD11b_gene_PU1 = control_parameter_dictionary["W_gene_CD11b_gene_PU1"]
	W_gene_CD14_RNAP = control_parameter_dictionary["W_gene_CD14_RNAP"]
	W_gene_CD14_gene_EGR1 = control_parameter_dictionary["W_gene_CD14_gene_EGR1"]
	W_gene_CD14_gene_CEBPa = control_parameter_dictionary["W_gene_CD14_gene_CEBPa"]
	W_gene_CD14_gene_PPARg = control_parameter_dictionary["W_gene_CD14_gene_PPARg"]
	W_gene_CD38_RNAP = control_parameter_dictionary["W_gene_CD38_RNAP"]
	W_gene_CD38_gene_RARa = control_parameter_dictionary["W_gene_CD38_gene_RARa"]
	W_gene_CD38_gene_PPARg = control_parameter_dictionary["W_gene_CD38_gene_PPARg"]
	W_gene_CD38_gene_IRF1 = control_parameter_dictionary["W_gene_CD38_gene_IRF1"]
	W_gene_CEBPa_RNAP = control_parameter_dictionary["W_gene_CEBPa_RNAP"]
	W_gene_CEBPa_gene_RARa = control_parameter_dictionary["W_gene_CEBPa_gene_RARa"]
	W_gene_CEBPa_gene_PPARg = control_parameter_dictionary["W_gene_CEBPa_gene_PPARg"]
	W_gene_CEBPa_gene_CEBPa = control_parameter_dictionary["W_gene_CEBPa_gene_CEBPa"]
	W_gene_CEBPa_gene_GFI1 = control_parameter_dictionary["W_gene_CEBPa_gene_GFI1"]
	W_gene_E2F_RNAP = control_parameter_dictionary["W_gene_E2F_RNAP"]
	W_gene_E2F_gene_E2F = control_parameter_dictionary["W_gene_E2F_gene_E2F"]
	W_gene_E2F_gene_PPARg = control_parameter_dictionary["W_gene_E2F_gene_PPARg"]
	W_gene_E2F_gene_CEBPa = control_parameter_dictionary["W_gene_E2F_gene_CEBPa"]
	W_gene_E2F_gene_GFI1 = control_parameter_dictionary["W_gene_E2F_gene_GFI1"]
	W_gene_EGR1_RNAP = control_parameter_dictionary["W_gene_EGR1_RNAP"]
	W_gene_EGR1_gene_RARa = control_parameter_dictionary["W_gene_EGR1_gene_RARa"]
	W_gene_EGR1_gene_PU1 = control_parameter_dictionary["W_gene_EGR1_gene_PU1"]
	W_gene_EGR1_gene_PPARg = control_parameter_dictionary["W_gene_EGR1_gene_PPARg"]
	W_gene_EGR1_gene_GFI1 = control_parameter_dictionary["W_gene_EGR1_gene_GFI1"]
	W_gene_GFI1_RNAP = control_parameter_dictionary["W_gene_GFI1_RNAP"]
	W_gene_GFI1_gene_CEBPa = control_parameter_dictionary["W_gene_GFI1_gene_CEBPa"]
	W_gene_GFI1_gene_EGR1 = control_parameter_dictionary["W_gene_GFI1_gene_EGR1"]
	W_gene_IRF1_RNAP = control_parameter_dictionary["W_gene_IRF1_RNAP"]
	W_gene_IRF1_gene_RARa = control_parameter_dictionary["W_gene_IRF1_gene_RARa"]
	W_gene_IRF1_gene_AhR = control_parameter_dictionary["W_gene_IRF1_gene_AhR"]
	W_gene_IRF1_gene_PPARg = control_parameter_dictionary["W_gene_IRF1_gene_PPARg"]
	W_gene_OCT1_RNAP = control_parameter_dictionary["W_gene_OCT1_RNAP"]
	W_gene_OCT1_gene_PPARg = control_parameter_dictionary["W_gene_OCT1_gene_PPARg"]
	W_gene_OCT4_RNAP = control_parameter_dictionary["W_gene_OCT4_RNAP"]
	W_gene_OCT4_gene_RARa = control_parameter_dictionary["W_gene_OCT4_gene_RARa"]
	W_gene_OCT4_gene_AhR = control_parameter_dictionary["W_gene_OCT4_gene_AhR"]
	W_gene_P21_RNAP = control_parameter_dictionary["W_gene_P21_RNAP"]
	W_gene_P21_gene_RARa = control_parameter_dictionary["W_gene_P21_gene_RARa"]
	W_gene_P21_gene_PPARg = control_parameter_dictionary["W_gene_P21_gene_PPARg"]
	W_gene_P21_gene_PU1 = control_parameter_dictionary["W_gene_P21_gene_PU1"]
	W_gene_P21_gene_IRF1 = control_parameter_dictionary["W_gene_P21_gene_IRF1"]
	W_gene_P21_gene_CEBPa = control_parameter_dictionary["W_gene_P21_gene_CEBPa"]
	W_gene_P21_gene_AP1 = control_parameter_dictionary["W_gene_P21_gene_AP1"]
	W_gene_P21_gene_GFI1 = control_parameter_dictionary["W_gene_P21_gene_GFI1"]
	W_gene_P47Phox_RNAP = control_parameter_dictionary["W_gene_P47Phox_RNAP"]
	W_gene_P47Phox_gene_PU1 = control_parameter_dictionary["W_gene_P47Phox_gene_PU1"]
	W_gene_P47Phox_gene_CEBPa = control_parameter_dictionary["W_gene_P47Phox_gene_CEBPa"]
	W_gene_P47Phox_gene_PPARg = control_parameter_dictionary["W_gene_P47Phox_gene_PPARg"]
	W_gene_PPARg_RNAP = control_parameter_dictionary["W_gene_PPARg_RNAP"]
	W_gene_PPARg_gene_CEBPa = control_parameter_dictionary["W_gene_PPARg_gene_CEBPa"]
	W_gene_PPARg_gene_EGR1 = control_parameter_dictionary["W_gene_PPARg_gene_EGR1"]
	W_gene_PPARg_gene_PU1 = control_parameter_dictionary["W_gene_PPARg_gene_PU1"]
	W_gene_PPARg_gene_AP1 = control_parameter_dictionary["W_gene_PPARg_gene_AP1"]
	W_gene_PU1_RNAP = control_parameter_dictionary["W_gene_PU1_RNAP"]
	W_gene_PU1_gene_RARa = control_parameter_dictionary["W_gene_PU1_gene_RARa"]
	W_gene_PU1_gene_AP1 = control_parameter_dictionary["W_gene_PU1_gene_AP1"]
	W_gene_PU1_gene_OCT1 = control_parameter_dictionary["W_gene_PU1_gene_OCT1"]
	W_gene_PU1_gene_AhR = control_parameter_dictionary["W_gene_PU1_gene_AhR"]
	W_gene_PU1_gene_GFI1 = control_parameter_dictionary["W_gene_PU1_gene_GFI1"]
	W_gene_RARa_RNAP = control_parameter_dictionary["W_gene_RARa_RNAP"]
	W_gene_RARa_gene_RARa = control_parameter_dictionary["W_gene_RARa_gene_RARa"]

	# Control function for gene_AP1 - 
	b_gene_AP1_gene_AhR = ((protein_gene_AhR)^(n_gene_AP1_gene_AhR))/(K_gene_AP1_gene_AhR^(n_gene_AP1_gene_AhR)+protein_gene_AhR^(n_gene_AP1_gene_AhR))
	b_gene_AP1_gene_PU1 = ((protein_gene_PU1)^(n_gene_AP1_gene_PU1))/(K_gene_AP1_gene_PU1^(n_gene_AP1_gene_PU1)+protein_gene_PU1^(n_gene_AP1_gene_PU1))
	b_gene_AP1_gene_PPARg = ((protein_gene_PPARg)^(n_gene_AP1_gene_PPARg))/(K_gene_AP1_gene_PPARg^(n_gene_AP1_gene_PPARg)+protein_gene_PPARg^(n_gene_AP1_gene_PPARg))
	control_array[1] = (W_gene_AP1_gene_AhR*b_gene_AP1_gene_AhR+W_gene_AP1_gene_PU1*b_gene_AP1_gene_PU1)/(1+W_gene_AP1_RNAP+W_gene_AP1_gene_AhR*b_gene_AP1_gene_AhR+W_gene_AP1_gene_PU1*b_gene_AP1_gene_PU1+W_gene_AP1_gene_PPARg*b_gene_AP1_gene_PPARg)

	# Control function for gene_AhR - 
	b_gene_AhR_gene_RARa = ((protein_gene_RARa)^(n_gene_AhR_gene_RARa))/(K_gene_AhR_gene_RARa^(n_gene_AhR_gene_RARa)+protein_gene_RARa^(n_gene_AhR_gene_RARa))
	control_array[2] = (W_gene_AhR_gene_RARa*b_gene_AhR_gene_RARa)/(1+W_gene_AhR_RNAP+W_gene_AhR_gene_RARa*b_gene_AhR_gene_RARa)

	# Control function for gene_CD11b - 
	b_gene_CD11b_gene_PU1 = ((protein_gene_PU1)^(n_gene_CD11b_gene_PU1))/(K_gene_CD11b_gene_PU1^(n_gene_CD11b_gene_PU1)+protein_gene_PU1^(n_gene_CD11b_gene_PU1))
	control_array[3] = (W_gene_CD11b_gene_PU1*b_gene_CD11b_gene_PU1)/(1+W_gene_CD11b_RNAP+W_gene_CD11b_gene_PU1*b_gene_CD11b_gene_PU1)

	# Control function for gene_CD14 - 
	b_gene_CD14_gene_EGR1 = ((protein_gene_EGR1)^(n_gene_CD14_gene_EGR1))/(K_gene_CD14_gene_EGR1^(n_gene_CD14_gene_EGR1)+protein_gene_EGR1^(n_gene_CD14_gene_EGR1))
	b_gene_CD14_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_CD14_gene_CEBPa))/(K_gene_CD14_gene_CEBPa^(n_gene_CD14_gene_CEBPa)+protein_gene_CEBPa^(n_gene_CD14_gene_CEBPa))
	b_gene_CD14_gene_PPARg = ((protein_gene_PPARg)^(n_gene_CD14_gene_PPARg))/(K_gene_CD14_gene_PPARg^(n_gene_CD14_gene_PPARg)+protein_gene_PPARg^(n_gene_CD14_gene_PPARg))
	control_array[4] = (W_gene_CD14_gene_EGR1*b_gene_CD14_gene_EGR1+W_gene_CD14_gene_CEBPa*b_gene_CD14_gene_CEBPa+W_gene_CD14_gene_PPARg*b_gene_CD14_gene_PPARg)/(1+W_gene_CD14_RNAP+W_gene_CD14_gene_EGR1*b_gene_CD14_gene_EGR1+W_gene_CD14_gene_CEBPa*b_gene_CD14_gene_CEBPa+W_gene_CD14_gene_PPARg*b_gene_CD14_gene_PPARg)

	# Control function for gene_CD38 - 
	b_gene_CD38_gene_RARa = ((protein_gene_RARa)^(n_gene_CD38_gene_RARa))/(K_gene_CD38_gene_RARa^(n_gene_CD38_gene_RARa)+protein_gene_RARa^(n_gene_CD38_gene_RARa))
	b_gene_CD38_gene_PPARg = ((protein_gene_PPARg)^(n_gene_CD38_gene_PPARg))/(K_gene_CD38_gene_PPARg^(n_gene_CD38_gene_PPARg)+protein_gene_PPARg^(n_gene_CD38_gene_PPARg))
	b_gene_CD38_gene_IRF1 = ((protein_gene_IRF1)^(n_gene_CD38_gene_IRF1))/(K_gene_CD38_gene_IRF1^(n_gene_CD38_gene_IRF1)+protein_gene_IRF1^(n_gene_CD38_gene_IRF1))
	control_array[5] = (W_gene_CD38_gene_RARa*b_gene_CD38_gene_RARa+W_gene_CD38_gene_PPARg*b_gene_CD38_gene_PPARg+W_gene_CD38_gene_IRF1*b_gene_CD38_gene_IRF1)/(1+W_gene_CD38_RNAP+W_gene_CD38_gene_RARa*b_gene_CD38_gene_RARa+W_gene_CD38_gene_PPARg*b_gene_CD38_gene_PPARg+W_gene_CD38_gene_IRF1*b_gene_CD38_gene_IRF1)

	# Control function for gene_CEBPa - 
	b_gene_CEBPa_gene_RARa = ((protein_gene_RARa)^(n_gene_CEBPa_gene_RARa))/(K_gene_CEBPa_gene_RARa^(n_gene_CEBPa_gene_RARa)+protein_gene_RARa^(n_gene_CEBPa_gene_RARa))
	b_gene_CEBPa_gene_PPARg = ((protein_gene_PPARg)^(n_gene_CEBPa_gene_PPARg))/(K_gene_CEBPa_gene_PPARg^(n_gene_CEBPa_gene_PPARg)+protein_gene_PPARg^(n_gene_CEBPa_gene_PPARg))
	b_gene_CEBPa_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_CEBPa_gene_CEBPa))/(K_gene_CEBPa_gene_CEBPa^(n_gene_CEBPa_gene_CEBPa)+protein_gene_CEBPa^(n_gene_CEBPa_gene_CEBPa))
	b_gene_CEBPa_gene_GFI1 = ((protein_gene_GFI1)^(n_gene_CEBPa_gene_GFI1))/(K_gene_CEBPa_gene_GFI1^(n_gene_CEBPa_gene_GFI1)+protein_gene_GFI1^(n_gene_CEBPa_gene_GFI1))
	control_array[6] = (W_gene_CEBPa_gene_RARa*b_gene_CEBPa_gene_RARa+W_gene_CEBPa_gene_PPARg*b_gene_CEBPa_gene_PPARg+W_gene_CEBPa_gene_CEBPa*b_gene_CEBPa_gene_CEBPa)/(1+W_gene_CEBPa_RNAP+W_gene_CEBPa_gene_RARa*b_gene_CEBPa_gene_RARa+W_gene_CEBPa_gene_PPARg*b_gene_CEBPa_gene_PPARg+W_gene_CEBPa_gene_CEBPa*b_gene_CEBPa_gene_CEBPa+W_gene_CEBPa_gene_GFI1*b_gene_CEBPa_gene_GFI1)

	# Control function for gene_E2F - 
	b_gene_E2F_gene_E2F = ((protein_gene_E2F)^(n_gene_E2F_gene_E2F))/(K_gene_E2F_gene_E2F^(n_gene_E2F_gene_E2F)+protein_gene_E2F^(n_gene_E2F_gene_E2F))
	b_gene_E2F_gene_PPARg = ((protein_gene_PPARg)^(n_gene_E2F_gene_PPARg))/(K_gene_E2F_gene_PPARg^(n_gene_E2F_gene_PPARg)+protein_gene_PPARg^(n_gene_E2F_gene_PPARg))
	b_gene_E2F_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_E2F_gene_CEBPa))/(K_gene_E2F_gene_CEBPa^(n_gene_E2F_gene_CEBPa)+protein_gene_CEBPa^(n_gene_E2F_gene_CEBPa))
	b_gene_E2F_gene_GFI1 = ((protein_gene_GFI1)^(n_gene_E2F_gene_GFI1))/(K_gene_E2F_gene_GFI1^(n_gene_E2F_gene_GFI1)+protein_gene_GFI1^(n_gene_E2F_gene_GFI1))
	control_array[7] = (W_gene_E2F_gene_E2F*b_gene_E2F_gene_E2F)/(1+W_gene_E2F_RNAP+W_gene_E2F_gene_E2F*b_gene_E2F_gene_E2F+W_gene_E2F_gene_PPARg*b_gene_E2F_gene_PPARg+W_gene_E2F_gene_CEBPa*b_gene_E2F_gene_CEBPa+W_gene_E2F_gene_GFI1*b_gene_E2F_gene_GFI1)

	# Control function for gene_EGR1 - 
	b_gene_EGR1_gene_RARa = ((protein_gene_RARa)^(n_gene_EGR1_gene_RARa))/(K_gene_EGR1_gene_RARa^(n_gene_EGR1_gene_RARa)+protein_gene_RARa^(n_gene_EGR1_gene_RARa))
	b_gene_EGR1_gene_PU1 = ((protein_gene_PU1)^(n_gene_EGR1_gene_PU1))/(K_gene_EGR1_gene_PU1^(n_gene_EGR1_gene_PU1)+protein_gene_PU1^(n_gene_EGR1_gene_PU1))
	b_gene_EGR1_gene_PPARg = ((protein_gene_PPARg)^(n_gene_EGR1_gene_PPARg))/(K_gene_EGR1_gene_PPARg^(n_gene_EGR1_gene_PPARg)+protein_gene_PPARg^(n_gene_EGR1_gene_PPARg))
	b_gene_EGR1_gene_GFI1 = ((protein_gene_GFI1)^(n_gene_EGR1_gene_GFI1))/(K_gene_EGR1_gene_GFI1^(n_gene_EGR1_gene_GFI1)+protein_gene_GFI1^(n_gene_EGR1_gene_GFI1))
	control_array[8] = (W_gene_EGR1_gene_RARa*b_gene_EGR1_gene_RARa+W_gene_EGR1_gene_PU1*b_gene_EGR1_gene_PU1)/(1+W_gene_EGR1_RNAP+W_gene_EGR1_gene_RARa*b_gene_EGR1_gene_RARa+W_gene_EGR1_gene_PU1*b_gene_EGR1_gene_PU1+W_gene_EGR1_gene_PPARg*b_gene_EGR1_gene_PPARg+W_gene_EGR1_gene_GFI1*b_gene_EGR1_gene_GFI1)

	# Control function for gene_GFI1 - 
	b_gene_GFI1_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_GFI1_gene_CEBPa))/(K_gene_GFI1_gene_CEBPa^(n_gene_GFI1_gene_CEBPa)+protein_gene_CEBPa^(n_gene_GFI1_gene_CEBPa))
	b_gene_GFI1_gene_EGR1 = ((protein_gene_EGR1)^(n_gene_GFI1_gene_EGR1))/(K_gene_GFI1_gene_EGR1^(n_gene_GFI1_gene_EGR1)+protein_gene_EGR1^(n_gene_GFI1_gene_EGR1))
	control_array[9] = (W_gene_GFI1_gene_CEBPa*b_gene_GFI1_gene_CEBPa)/(1+W_gene_GFI1_RNAP+W_gene_GFI1_gene_CEBPa*b_gene_GFI1_gene_CEBPa+W_gene_GFI1_gene_EGR1*b_gene_GFI1_gene_EGR1)

	# Control function for gene_IRF1 - 
	b_gene_IRF1_gene_RARa = ((protein_gene_RARa)^(n_gene_IRF1_gene_RARa))/(K_gene_IRF1_gene_RARa^(n_gene_IRF1_gene_RARa)+protein_gene_RARa^(n_gene_IRF1_gene_RARa))
	b_gene_IRF1_gene_AhR = ((protein_gene_AhR)^(n_gene_IRF1_gene_AhR))/(K_gene_IRF1_gene_AhR^(n_gene_IRF1_gene_AhR)+protein_gene_AhR^(n_gene_IRF1_gene_AhR))
	b_gene_IRF1_gene_PPARg = ((protein_gene_PPARg)^(n_gene_IRF1_gene_PPARg))/(K_gene_IRF1_gene_PPARg^(n_gene_IRF1_gene_PPARg)+protein_gene_PPARg^(n_gene_IRF1_gene_PPARg))
	control_array[10] = (W_gene_IRF1_gene_RARa*b_gene_IRF1_gene_RARa+W_gene_IRF1_gene_AhR*b_gene_IRF1_gene_AhR+W_gene_IRF1_gene_PPARg*b_gene_IRF1_gene_PPARg)/(1+W_gene_IRF1_RNAP+W_gene_IRF1_gene_RARa*b_gene_IRF1_gene_RARa+W_gene_IRF1_gene_AhR*b_gene_IRF1_gene_AhR+W_gene_IRF1_gene_PPARg*b_gene_IRF1_gene_PPARg)

	# Control function for gene_OCT1 - 
	b_gene_OCT1_gene_PPARg = ((protein_gene_PPARg)^(n_gene_OCT1_gene_PPARg))/(K_gene_OCT1_gene_PPARg^(n_gene_OCT1_gene_PPARg)+protein_gene_PPARg^(n_gene_OCT1_gene_PPARg))
	control_array[11] = (W_gene_OCT1_gene_PPARg*b_gene_OCT1_gene_PPARg)/(1+W_gene_OCT1_RNAP+W_gene_OCT1_gene_PPARg*b_gene_OCT1_gene_PPARg)

	# Control function for gene_OCT4 - 
	b_gene_OCT4_gene_RARa = ((protein_gene_RARa)^(n_gene_OCT4_gene_RARa))/(K_gene_OCT4_gene_RARa^(n_gene_OCT4_gene_RARa)+protein_gene_RARa^(n_gene_OCT4_gene_RARa))
	b_gene_OCT4_gene_AhR = ((protein_gene_AhR)^(n_gene_OCT4_gene_AhR))/(K_gene_OCT4_gene_AhR^(n_gene_OCT4_gene_AhR)+protein_gene_AhR^(n_gene_OCT4_gene_AhR))
	control_array[12] = ()/(1+W_gene_OCT4_RNAP+W_gene_OCT4_gene_RARa*b_gene_OCT4_gene_RARa+W_gene_OCT4_gene_AhR*b_gene_OCT4_gene_AhR)

	# Control function for gene_P21 - 
	b_gene_P21_gene_RARa = ((protein_gene_RARa)^(n_gene_P21_gene_RARa))/(K_gene_P21_gene_RARa^(n_gene_P21_gene_RARa)+protein_gene_RARa^(n_gene_P21_gene_RARa))
	b_gene_P21_gene_PPARg = ((protein_gene_PPARg)^(n_gene_P21_gene_PPARg))/(K_gene_P21_gene_PPARg^(n_gene_P21_gene_PPARg)+protein_gene_PPARg^(n_gene_P21_gene_PPARg))
	b_gene_P21_gene_PU1 = ((protein_gene_PU1)^(n_gene_P21_gene_PU1))/(K_gene_P21_gene_PU1^(n_gene_P21_gene_PU1)+protein_gene_PU1^(n_gene_P21_gene_PU1))
	b_gene_P21_gene_IRF1 = ((protein_gene_IRF1)^(n_gene_P21_gene_IRF1))/(K_gene_P21_gene_IRF1^(n_gene_P21_gene_IRF1)+protein_gene_IRF1^(n_gene_P21_gene_IRF1))
	b_gene_P21_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_P21_gene_CEBPa))/(K_gene_P21_gene_CEBPa^(n_gene_P21_gene_CEBPa)+protein_gene_CEBPa^(n_gene_P21_gene_CEBPa))
	b_gene_P21_gene_AP1 = ((protein_gene_AP1)^(n_gene_P21_gene_AP1))/(K_gene_P21_gene_AP1^(n_gene_P21_gene_AP1)+protein_gene_AP1^(n_gene_P21_gene_AP1))
	b_gene_P21_gene_GFI1 = ((protein_gene_GFI1)^(n_gene_P21_gene_GFI1))/(K_gene_P21_gene_GFI1^(n_gene_P21_gene_GFI1)+protein_gene_GFI1^(n_gene_P21_gene_GFI1))
	control_array[13] = (W_gene_P21_gene_RARa*b_gene_P21_gene_RARa+W_gene_P21_gene_PPARg*b_gene_P21_gene_PPARg+W_gene_P21_gene_PU1*b_gene_P21_gene_PU1+W_gene_P21_gene_IRF1*b_gene_P21_gene_IRF1+W_gene_P21_gene_CEBPa*b_gene_P21_gene_CEBPa+W_gene_P21_gene_AP1*b_gene_P21_gene_AP1)/(1+W_gene_P21_RNAP+W_gene_P21_gene_RARa*b_gene_P21_gene_RARa+W_gene_P21_gene_PPARg*b_gene_P21_gene_PPARg+W_gene_P21_gene_PU1*b_gene_P21_gene_PU1+W_gene_P21_gene_IRF1*b_gene_P21_gene_IRF1+W_gene_P21_gene_CEBPa*b_gene_P21_gene_CEBPa+W_gene_P21_gene_AP1*b_gene_P21_gene_AP1+W_gene_P21_gene_GFI1*b_gene_P21_gene_GFI1)

	# Control function for gene_P47Phox - 
	b_gene_P47Phox_gene_PU1 = ((protein_gene_PU1)^(n_gene_P47Phox_gene_PU1))/(K_gene_P47Phox_gene_PU1^(n_gene_P47Phox_gene_PU1)+protein_gene_PU1^(n_gene_P47Phox_gene_PU1))
	b_gene_P47Phox_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_P47Phox_gene_CEBPa))/(K_gene_P47Phox_gene_CEBPa^(n_gene_P47Phox_gene_CEBPa)+protein_gene_CEBPa^(n_gene_P47Phox_gene_CEBPa))
	b_gene_P47Phox_gene_PPARg = ((protein_gene_PPARg)^(n_gene_P47Phox_gene_PPARg))/(K_gene_P47Phox_gene_PPARg^(n_gene_P47Phox_gene_PPARg)+protein_gene_PPARg^(n_gene_P47Phox_gene_PPARg))
	control_array[14] = (W_gene_P47Phox_gene_PU1*b_gene_P47Phox_gene_PU1+W_gene_P47Phox_gene_CEBPa*b_gene_P47Phox_gene_CEBPa)/(1+W_gene_P47Phox_RNAP+W_gene_P47Phox_gene_PU1*b_gene_P47Phox_gene_PU1+W_gene_P47Phox_gene_CEBPa*b_gene_P47Phox_gene_CEBPa+W_gene_P47Phox_gene_PPARg*b_gene_P47Phox_gene_PPARg)

	# Control function for gene_PPARg - 
	b_gene_PPARg_gene_CEBPa = ((protein_gene_CEBPa)^(n_gene_PPARg_gene_CEBPa))/(K_gene_PPARg_gene_CEBPa^(n_gene_PPARg_gene_CEBPa)+protein_gene_CEBPa^(n_gene_PPARg_gene_CEBPa))
	b_gene_PPARg_gene_EGR1 = ((protein_gene_EGR1)^(n_gene_PPARg_gene_EGR1))/(K_gene_PPARg_gene_EGR1^(n_gene_PPARg_gene_EGR1)+protein_gene_EGR1^(n_gene_PPARg_gene_EGR1))
	b_gene_PPARg_gene_PU1 = ((protein_gene_PU1)^(n_gene_PPARg_gene_PU1))/(K_gene_PPARg_gene_PU1^(n_gene_PPARg_gene_PU1)+protein_gene_PU1^(n_gene_PPARg_gene_PU1))
	b_gene_PPARg_gene_AP1 = ((protein_gene_AP1)^(n_gene_PPARg_gene_AP1))/(K_gene_PPARg_gene_AP1^(n_gene_PPARg_gene_AP1)+protein_gene_AP1^(n_gene_PPARg_gene_AP1))
	control_array[15] = (W_gene_PPARg_gene_CEBPa*b_gene_PPARg_gene_CEBPa+W_gene_PPARg_gene_EGR1*b_gene_PPARg_gene_EGR1)/(1+W_gene_PPARg_RNAP+W_gene_PPARg_gene_CEBPa*b_gene_PPARg_gene_CEBPa+W_gene_PPARg_gene_EGR1*b_gene_PPARg_gene_EGR1+W_gene_PPARg_gene_PU1*b_gene_PPARg_gene_PU1+W_gene_PPARg_gene_AP1*b_gene_PPARg_gene_AP1)

	# Control function for gene_PU1 - 
	b_gene_PU1_gene_RARa = ((protein_gene_RARa)^(n_gene_PU1_gene_RARa))/(K_gene_PU1_gene_RARa^(n_gene_PU1_gene_RARa)+protein_gene_RARa^(n_gene_PU1_gene_RARa))
	b_gene_PU1_gene_AP1 = ((protein_gene_AP1)^(n_gene_PU1_gene_AP1))/(K_gene_PU1_gene_AP1^(n_gene_PU1_gene_AP1)+protein_gene_AP1^(n_gene_PU1_gene_AP1))
	b_gene_PU1_gene_OCT1 = ((protein_gene_OCT1)^(n_gene_PU1_gene_OCT1))/(K_gene_PU1_gene_OCT1^(n_gene_PU1_gene_OCT1)+protein_gene_OCT1^(n_gene_PU1_gene_OCT1))
	b_gene_PU1_gene_AhR = ((protein_gene_AhR)^(n_gene_PU1_gene_AhR))/(K_gene_PU1_gene_AhR^(n_gene_PU1_gene_AhR)+protein_gene_AhR^(n_gene_PU1_gene_AhR))
	b_gene_PU1_gene_GFI1 = ((protein_gene_GFI1)^(n_gene_PU1_gene_GFI1))/(K_gene_PU1_gene_GFI1^(n_gene_PU1_gene_GFI1)+protein_gene_GFI1^(n_gene_PU1_gene_GFI1))
	control_array[16] = (W_gene_PU1_gene_RARa*b_gene_PU1_gene_RARa+W_gene_PU1_gene_AP1*b_gene_PU1_gene_AP1+W_gene_PU1_gene_OCT1*b_gene_PU1_gene_OCT1)/(1+W_gene_PU1_RNAP+W_gene_PU1_gene_RARa*b_gene_PU1_gene_RARa+W_gene_PU1_gene_AP1*b_gene_PU1_gene_AP1+W_gene_PU1_gene_OCT1*b_gene_PU1_gene_OCT1+W_gene_PU1_gene_AhR*b_gene_PU1_gene_AhR+W_gene_PU1_gene_GFI1*b_gene_PU1_gene_GFI1)

	# Control function for gene_RARa - 
	b_gene_RARa_gene_RARa = ((protein_gene_RARa)^(n_gene_RARa_gene_RARa))/(K_gene_RARa_gene_RARa^(n_gene_RARa_gene_RARa)+protein_gene_RARa^(n_gene_RARa_gene_RARa))
	control_array[17] = (W_gene_RARa_gene_RARa*b_gene_RARa_gene_RARa)/(1+W_gene_RARa_RNAP+W_gene_RARa_gene_RARa*b_gene_RARa_gene_RARa)


	return control_array
end
