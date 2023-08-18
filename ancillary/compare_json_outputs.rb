#!/usr/bin/env ruby

require 'json'

content_ref = JSON.parse(File.read(ARGV[0]))
content_new = JSON.parse(File.read(ARGV[1]))

content_ref["cross sections"].each do |key, _|
  content_ref["cross sections"][key].each do |key2, _|
    xsec_ref = content_ref["cross sections"][key][key2]["res"]
    err_ref  = content_ref["cross sections"][key][key2]["err"]
    xsec_new = content_new["cross sections"][key][key2]["res"]
    err_new  = content_new["cross sections"][key][key2]["err"]

    exit 1 if err_ref < 0 || err_new < 0

    if ([err_ref, err_new].max > 0)
      sigma = (xsec_ref-xsec_new).abs/Math.sqrt(err_ref**2 + err_new**2)
      puts "#{key} #{key2} difference:".ljust(40, ' ') + "#{sprintf('%.1f', sigma)}σ      ".rjust(12, ' ') + "(#{xsec_ref} +/- #{err_ref} vs. #{xsec_new} +/- #{err_new})"

      exit 1 if sigma > 2.0
    end
  end
end
