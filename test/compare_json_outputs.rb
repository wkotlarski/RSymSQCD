#!/usr/bin/env ruby

require 'json'
require 'securerandom'

uuid = SecureRandom.uuid
system({"CUBACORES" => ARGV[3]}, "#{ARGV[0]}/RSymSQCD --card #{ARGV[1]} --json-outputfile-name=#{uuid}")
exit 1 if $?.exitstatus != 0

file_ref     = File.read(ARGV[2])
file_to_test = File.read(uuid)
content_ref = JSON.parse(file_ref)
content_new = JSON.parse(file_to_test)

content_ref["cross sections"].each do |key, _|
  content_ref["cross sections"][key].each do |key2, _|
    xsec_ref = content_ref["cross sections"][key][key2]["res"]
    err_ref  = content_ref["cross sections"][key][key2]["err"]
    xsec_new = content_new["cross sections"][key][key2]["res"]
    err_new  = content_new["cross sections"][key][key2]["err"]

    exit 1 if err_ref < 0 || err_new < 0

    sigma = (xsec_ref-xsec_new).abs/Math.sqrt(err_ref**2 + err_new**2)
    puts "#{key} #{key2}: differ by #{sprintf('%.1f', sigma)}σ (#{xsec_ref} +/- #{err_ref} vs. #{xsec_new} +/- #{err_new})"

    exit 1 if sigma > 2.0
  end
end
