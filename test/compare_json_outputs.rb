#!/usr/bin/env ruby

require 'json'
require 'securerandom'

uuid = SecureRandom.uuid
system("CUBACORES=3 #{ARGV[0]}/RSymSQCD --card #{ARGV[1]} --json-outputfile-name=#{uuid}")

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

    puts "#{key} #{key2}: #{xsec_ref} v. #{xsec_new} (1-ratio: #{1-xsec_ref/xsec_new}, max_error: #{[err_ref, err_new].max}"

    exit 1 if (xsec_ref-xsec_new).abs > 3.0*Math.sqrt(err_ref**2 + err_new**2)
  end
end
