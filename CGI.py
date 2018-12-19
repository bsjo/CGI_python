#!/usr/bin/env python

import requests, StringIO, zipfile, os, datetime, re, string


class CGI_parsing():
	def __init__(self):
		print("START\n")
		temp_today = datetime.datetime.today()
		self.today = str(temp_today).split(" ")[0]
		print(temp_today)

	def CGI_download(self, url, input_path):
		data = requests.get(url)
		f = StringIO.StringIO()
		f.write(data.content)
		input_zip = zipfile.ZipFile(f)
		r_data = [input_zip.read(i) for i in input_zip.namelist() if i == 'cgi_biomarkers_per_variant.tsv'][0]
		self.r_list = r_data.splitlines()
		#CGI making
		if not os.path.exists(input_path): os.makedirs(input_path)
		CGI_down = '{0}/CGI_CancerBiomarkers_{1}'.format(input_path, self.today)
		with open(CGI_down, 'w') as w:
			w.write(r_data)

	def CGI_processing(self, output_path):
		if not os.path.exists(output_path): os.makedirs(output_path)
		output_path_n = os.path.join(output_path, "Parsed_CGI_{0}".format(self.today))
		with open(output_path_n, "w") as w:
			header = '\t'.join(['Chromosome', 'Position', 'G_change', 'P_change', 'Exon_order', 'Gene', 'Transcript', 'Disease', 'Drug', 'Drug_status', 'Evidence']) + '\n'
			w.write(header)
			for line in self.r_list[1:]:
				line_split = line.split('\t')
				Genomic_pos = line_split[21]
				Drug = line_split[8] ; Drug_status = line_split[11] ; Drug_evi = line_split[12]
				Cancer = line_split[16] ; Gene = line_split[13] ; Transcript = line_split[26]
				out_result = ""
				if len(Genomic_pos) > 0:
					chrom = re.findall(r'chr\w+', Genomic_pos)[0] ; exon_region = re.findall(r'exon_\d+', line_split[24])[0]
					G_pos = re.findall(r'g.\d+', Genomic_pos)[0].split(".")[1]
					G_change = "" ; P_change = line_split[22].split(":")[1]
					if re.search('del', Genomic_pos):
						G_temp = re.findall('[a-z]+[A-Z]*',  Genomic_pos)
						G_change  = G_temp[2].translate(None,string.ascii_lowercase) + ">" + G_temp[3].translate(None,string.ascii_lowercase)
					else:
						G_change = re.findall(r'\D+>\D+', Genomic_pos)[0]
					out_result = '\t'.join([chrom, G_pos, G_change, P_change, exon_region, Gene, Transcript, Cancer, Drug, Drug_status, Drug_evi]) + '\n'
				else: #No genomic pos
					chrom = G_pos = G_change = P_change = exon_region = '' ;
					out_result = '\t'.join([chrom, G_pos, G_change, P_change, exon_region, Gene, Transcript, Cancer, Drug, Drug_status, Drug_evi]) + '\n'
				w.write(out_result)

if __name__ == '__main__':

	url = "https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_latest.zip"
	input_path = 'Input'
	output_path = 'Output'

	#CGI_Class
	go = CGI_parsing()
	#CGI_CancerBiomarkers_Download_lastest
	go.CGI_download(url, input_path)
	#CGI_Parsing
	go.CGI_processing(output_path)

	print("\nEND")