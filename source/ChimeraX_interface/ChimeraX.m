function response = ChimeraX(command)

url = sprintf('http://127.0.0.1:51051/run?command=%s',urlencode(command));

options = weboptions('Timeout',60);

response = webread(url,options);