import subprocess
import re
from tqdm import tqdm

inv_pat = re.compile(r"recovery by inverse ok \? (\w+)")
trick_pat = re.compile(r"recovery by trick ok \? (\w+)")


def invoke(path, args=None):
    if args is None:
        args = []
    proc = subprocess.Popen([path] + args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.wait()
    res = proc.communicate()[0].decode("ascii")
    inv_match = inv_pat.findall(res)
    trick_match = trick_pat.findall(res)
    if not inv_match or not trick_match:
        print(res)
        print("nope")
        exit(1)
    return inv_match[0] == "true", trick_match[0] == "true"


def invoke_repeat(path, ntimes=1, args=None):
    inv_success, trick_success = 0, 0
    bar = tqdm(total=ntimes)
    for _ in range(ntimes):
        inv_ok, trick_ok = invoke(path, args)
        inv_success += inv_ok
        trick_success += trick_ok
        bar.update(1)
    bar.close()
    return inv_success, trick_success


if __name__ == "__main__":
    total = 100
    inv, trick = invoke_repeat("/home/msh/CKKS/cmake-build-release/heaanattack", ntimes=total)
    print("inv ok: {inv}/{total}, trick ok: {trick}/{total}".format(inv=inv, total=total, trick=trick))
