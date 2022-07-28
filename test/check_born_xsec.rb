require 'securerandom'
require 'json'

uuid = SecureRandom.uuid

system("#{ARGV[2]}/RSymSQCD --card #{ARGV[0]} --json-outputfile-name=#{uuid}")

file = File.read(uuid)
data = JSON.parse(file)

xsec = 0.0
err = 0.0
data["cross sections"].each do |e, v|
  xsec += v["tree"]["res"]
  err += (v["tree"]["err"])**2
end

exit 1 if (xsec-ARGV[1].to_f).abs > Math.sqrt(err)
