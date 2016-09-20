# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "hashicorp/precise64"
  config.vm.synced_folder ".", "/vagrant"
  config.vm.hostname = "cibox"
  config.vm.provision "shell",
    inline: "sudo apt-get update && sudo apt-get install make && cd /vagrant && make install-slurm"
end
