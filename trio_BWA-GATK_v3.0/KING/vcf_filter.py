import sys
from vcf import VCFReader, VCFWriter
from vcf.model import make_calldata_tuple


# usage: python vcf_filter.py test.vcf tmp.vcf

def make_hom_calldata(format, data):
    CallData =make_calldata_tuple(format.split(':'))
    call_data ={k: getattr(data, k) for k in format.split(':')}
    call_data.update(GT='0/0')
    return CallData(**call_data)

if __name__ == '__main__':
    vcf_in, vcf_out = sys.argv[1:]
    
    AF_thres = 0.30
    DP_thres = 20
    GQ_thres = 50
    QUAL_thres = 100

    vcf_reader =VCFReader(filename=vcf_in)
    with open(vcf_out, 'w') as fp:
        vcf_writer = VCFWriter(fp, vcf_reader)
        for record in vcf_reader:
            for sample in record.samples:
                if sample.data.AD and sample.data.DP:
                    af = sample.data.AD[1] / float(sample.data.DP)
                    if record.QUAL < QUAL_thres or sample.data.DP < DP_thres or sample.data.GQ < GQ_thres or af < AF_thres :
                        sample.data = make_hom_calldata(record.FORMAT, sample.data)
            vcf_writer.write_record(record)
