"""Read the ee fit handoff without extracting or executing archive members."""
import hashlib
import math
import re
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path


def read_archive(path):
    path = Path(path)
    with tarfile.open(path, 'r:gz') as archive:
        names = archive.getnames()
        def read(name):
            member = archive.getmember('Zee_fit/' + name)
            if not member.isfile() or member.size > 10_000_000:
                raise ValueError(f'Unsupported archive member: {member.name}')
            with archive.extractfile(member) as stream:
                return stream.read().decode('utf-8')
        match = re.search(r'^mu_signal\s+(\S+)\s+\+(\S+)\s+-(\S+)',
                          read('Fits/Zee_fit.txt'), re.MULTILINE)
        if not match:
            raise ValueError('Archived ee fit is missing mu_signal')
        value, up, down = map(float, match.groups())
        if not all(math.isfinite(v) for v in (value, up, down)) or min(up, down) < 0:
            raise ValueError('Invalid archived ee POI')
        errors = {}
        for line in read('Fits/Zee_fit_errDecomp_mu_signal.txt').splitlines():
            fields = line.split()
            if len(fields) == 4 and fields[0] in {'TOT_ERROR', 'STAT_ERROR', 'SYST_ERROR', 'MCSTAT_ERROR'}:
                errors[fields[0]] = dict(symmetric=float(fields[1]), up=float(fields[2]), down=float(fields[3]))
        measurement = ET.fromstring(read('RooStats/Zee_fit.xml')).find('Measurement')
        channel = ET.fromstring(read('RooStats/Zee_fit_SR.xml'))
        samples = [dict(name=s.attrib['Name'],
                        mc_stat=s.find('StatError') is not None and s.find('StatError').get('Activate') == 'True',
                        normfactors=[n.attrib for n in s.findall('NormFactor')])
                   for s in channel.findall('Sample')]
        systematic_nodes = [n.tag for n in channel.iter() if n.tag in {'OverallSys','HistoSys','ShapeSys','ShapeFactor'}]
        return dict(archive=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    poi='mu_signal', poi_value=value, poi_error_up=up, poi_error_down=down,
                    error_decomposition=errors, region=channel.get('Name'), samples=samples,
                    luminosity_model=measurement.attrib,
                    systematic_nodes=systematic_nodes,
                    workspace_available='Zee_fit/RooStats/Zee_fit_combined_Zee_fit_model.root' in names,
                    uncertainty_scope='Data and MC statistics only' if not systematic_nodes and measurement.get('LumiRelErr') == '0' else 'Inspect archived model',
                    use='Reference result only; rebuild from nominal datasets/z-ee histograms for common POI and external refits',
                    sigma_reference_60_120_pb=None)
