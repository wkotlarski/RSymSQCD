#!/usr/bin/env ruby

require 'json'

file_ref = File.read(ARGV[0])
file_new = File.read(ARGV[1])
data_hash_ref = JSON.parse(file_ref)
data_hash_new = JSON.parse(file_new)

data_hash_ref["cross sections"].each do |key, _|
  data_hash_ref["cross sections"][key].each do |key2, _|
    xsec_ref = data_hash_ref["cross sections"][key][key2]["res"]
    err_ref  = data_hash_ref["cross sections"][key][key2]["err"]
    xsec_new = data_hash_new["cross sections"][key][key2]["res"]
    err_new  = data_hash_new["cross sections"][key][key2]["err"]
    puts "Channel #{key} #{key2} difference: #{xsec_ref} +/- #{err_ref} vs. #{xsec_new} vs. #{err_new} (\#σ = #{(xsec_ref-xsec_new).abs/(3.0*Math.sqrt(err_ref**2 + err_new**2))})"
    if (xsec_ref-xsec_new).abs > 3.0*Math.sqrt(err_ref**2 + err_new**2)
      exit(1)
    end
  end
end
