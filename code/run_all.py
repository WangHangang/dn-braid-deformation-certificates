"""Run the complete audit from any working directory, Python >=3.10, no packages."""
from pathlib import Path
import argparse,subprocess,sys

root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--include-independent-audit',action='store_true',help='Also replay received symbolic/modular audit; needs SymPy.')
args=parser.parse_args()
commands=[['check_normal_form.py'],['check_uniform_identities.py','--ranks','4','5','6'],
          ['check_family.py'],['audit_certificates.py',str(root/'input_certificates')],
          ['check_prior_art_dictionary.py']]
if args.include_independent_audit:commands.append(['replay_independent_audit.py'])
for command in commands:
    print('RUN',*command,flush=True)
    subprocess.run([sys.executable,str(root/'code'/command[0]),*command[1:]],check=True,cwd=root)
print('ALL AUDITS PASSED',flush=True)
