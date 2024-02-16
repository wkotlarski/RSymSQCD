#!/usr/bin/env python3
import json
import sys
from math import sqrt

with open(sys.argv[1], "r") as read_file:
   data = json.load(read_file)

xsec = 0.
err = 0.
max_channel_error = 0.
max_pval = 0.
for subprocess_key, subprocess_value in data["cross sections"].items():
    for val in subprocess_value.items():
        xsec += val[1]['res']
        err  += (val[1]['err'])**2
        if val[1]['res'] > 0.:
            max_channel_error = max(max_channel_error, val[1]['err']/val[1]['res'])
            max_pval = max(max_pval, val[1]['p-val'])

err = sqrt(err)
print(f'Max channel error (%): {100.*max_channel_error}')
print(f'Max p-value: {max_pval}')

print(f'{data["masses"]["gluino"]}    {data["masses"]["squarks"]}    {data["mu_f"]}    {data["mu_r"]}    {xsec}    {err}    ({100.*err/xsec}%)')
