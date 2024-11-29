from bottle import Bottle, request, response, run, static_file, HTTPResponse, redirect
import os
#import uuid
import subprocess
import hashlib
from datetime import datetime

app = Bottle()

# Directory where job output folders will be stored
BASE_OUTPUT_DIR = '/storage/Data/AF2jobs'  # Change this to your desired path
log_file = os.path.join(BASE_OUTPUT_DIR, 'AF2.log')

@app.route('/favicon.ico', name='get_favicon')
def get_favicon():
    response.content_type = 'image/png'
    return "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAA+UlEQVQ4T6XSMQ5EQBQG4F8l0aiVCpE4gAIFcQaVwiX0GjeQOAAKpxARrVIpDqByANndeckmNsGO3Wllvnn/+wmP18EfR2CAqqoIwxBpmhKVJAmKosA8z19pAmzbhiiKaJqGLnieh23b0HUdH+A4DizLQhAEdKGuawzDgLZt+YEsy1CWJdhKWJw4ju8Bfd/D930anb0sSRJM06QJDMNAnueH09AOWAQGVFVFQBRFkGUZ67ryR2DA/uyBq5ZogqNndF2Hoij0aZomaJp22NIp8I7FgKuWuIGzlrgijOOIZVk+WnJdl9o6BfZ7OWrpJ2CP3gKufgauCFfAE09Z1dH0wbq/AAAAAElFTkSuQmCC"

@app.route('/AF2IBS/previous_results', name="previous_results")
def previous_results():
    log_lines = ''
    if os.path.isfile(log_file):
        with open(log_file, 'r') as f:
            log_line = f.read()
            if(log_line != None):
                log_lines += log_line
    pagestr = '<!DOCTYPE html>'
    pagestr += '<html>'
    pagestr += '<head>'
    pagestr += '<style type="text/css" media="screen">'
    pagestr += 'th { position:sticky; top:0px; background:white;}'
    pagestr += '</style>'
    pagestr += '</head>'
    pagestr += '<body>'
    pagestr += '<p><table style="width:100%">'
    pagestr += '<th>Date-time</th><th>User</th><th>IP</th><th>Job id</th><th>#tokens</th>'
    pagestr += log_lines
    pagestr += '</table></p>'
    pagestr += '</body>'
    pagestr += '</html>'
    return pagestr

@app.route('/AF2IBS', method=['GET', 'POST'])
def index():
    client_ip = request.environ.get('HTTP_X_FORWARDED_FOR') or request.environ.get('REMOTE_ADDR')
    log_lines = ''
    if os.path.isfile(log_file):
        with open(log_file, 'r') as f:
            log_line = f.read()
            if(log_line != None):
                log_lines += log_line
    if request.method == 'POST':
        # Get text parameters from the form
        sequences = request.forms.get('sequences').strip()+'\n'
        email = request.forms.get('email')
        tokens_var = request.forms.get('tokens_var')

        # Generate a unique job ID
        job_id = hashlib.md5(sequences.encode("utf-8")).hexdigest()

        # Create an output directory for the job
        job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)
        os.makedirs(job_output_dir, exist_ok=True)

        # Save sequences to a file inside job_output_dir
        sequences_file = os.path.join(job_output_dir, 'input_'+job_id+'.fasta')
        with open(sequences_file, 'w') as f:
            f.write(sequences)

        subprocess.call(['chmod', '664', sequences_file])
        # Launch the subprocess asynchronously
        my_env = os.environ.copy()
        my_env["IBSJOBNAME"] = job_id
        command = ['/storage/Alphafold/scripts/alphafold2_caller.bin', sequences_file, '3']
        subprocess.Popen(command, cwd=job_output_dir, env=my_env)

        # Return the job URL to the client
        job_url = request.urlparts.scheme + "://" + request.urlparts.netloc + app.get_url('job_results', job_id=job_id)
        with open(log_file, 'a') as f:
            f.write('<tr><td>'+str(datetime.now())+'</td><td>'+email.split("@")[0]+'</td><td>'+client_ip+'</td><td><a href="'+job_url+'">'+job_id+'</a></td><td>'+tokens_var+'</td></tr>\n')
        log_lines +='<tr><td>'+str(datetime.now())+'</td><td>'+email.split("@")[0]+'</td><td>'+client_ip+'</td><td><a href="'+job_url+'">'+job_id+'</a></td><td>'+tokens_var+'</td></tr>\n'
        return f'''
            <h1>Your Alphafold 2 job has been submitted</h1>
            <p>You can check the result at: <a href="{job_url}">{job_url}</a></p>
            <p>Keep the link above for future access.</p>
        '''

    # If GET request, display the submission form
    # Get the output of 'squeue'
    try:
        squeue_output = subprocess.check_output(['squeue'], universal_newlines=True)
    except subprocess.CalledProcessError as e:
        squeue_output = f"Failed to retrieve squeue output: {e}"
    pagestr = '''<!DOCTYPE html>
<html>
<head>
<script>
function countPipeInSubstring(str, start, end){
    let substr = str.slice(start, end);
    let regex = new RegExp(/[|]/, 'g');
    let matches = substr.match(regex);
    return matches ? matches.length : 0;
}
function checkValue() {
    if(sequences.value[0] != '>'){
        sequences.setCustomValidity("The first letter of the first line has to be the symbol >");
        return false;
    }
    if(sequences.value.indexOf('\\n') == -1){
        sequences.setCustomValidity("At least two lines. One for the sequence name, other for the sequence itself.");
        return false;
    }
    let str = sequences.value;
    let idx = str.indexOf('\\n');













    let yes_tokens = 0;
    let line = 1;
    let mymatch = null;
    while(idx != -1){
        str = str.slice(idx+1);
        idx = str.indexOf('\\n');
        if(str.slice(0,idx+1) != '' && str.slice(0,idx+1)[0] != '>'){
            yes_tokens += idx;
            let iJ = str.slice(0,idx+1).toUpperCase().indexOf('J');
            if(iJ != -1){
                sequences.setCustomValidity(`Invalid letter J at line ${line+1}, column ${iJ+1}`);
                return false;
            }
            let iO = str.slice(0,idx+1).toUpperCase().indexOf('O');
            if(iO != -1){
                sequences.setCustomValidity(`Invalid letter O at line ${line+1}, colunm ${iO+1}`);
                return false;
            }
            mymatch = str.slice(0,idx+1).match(/[0-9]/);
            if(mymatch){
                let iN = str.slice(0,idx+1).indexOf(mymatch[0]);
                sequences.setCustomValidity(`Invalid numeric value (${mymatch[0]}) at line ${line+1}, colunm ${iN+1}`);
                return false;
            }
            mymatch = str.slice(0,idx+1).match(/[^a-zA-Z||\\n\\r]/);
            if(mymatch){
                let iN = str.slice(0,idx+1).indexOf(mymatch[0]);
                sequences.setCustomValidity(`Invalid symbol (${mymatch[0]}) at line ${line+1}, colunm ${iN+1}`);
                return false;
            }
        }
        line++;
    }
    //last line exception
    if(str != '' && str[0] != '>'){
        yes_tokens += str.length;
        let iJ = str.toUpperCase().indexOf('J');
        if(iJ != -1){
            sequences.setCustomValidity(`Invalid letter J at line ${line}, column ${iJ+1}`);
            return false;
        }
        let iO = str.toUpperCase().indexOf('O');
        if(iO != -1){
            sequences.setCustomValidity(`Invalid letter O at line ${line}, colunm ${iO+1}`);
            return false;
        }
        mymatch = str.match(/[0-9]/);
        if(mymatch){
            let iN = str.indexOf(mymatch[0]);
            sequences.setCustomValidity(`Invalid numeric value (${mymatch[0]}) at line ${line}, colunm ${iN+1}`);
            return false;
        }
        mymatch = str.match(/[^a-zA-Z||\\n\\r]/);
        if(mymatch){
            let iN = str.indexOf(mymatch[0]);
            sequences.setCustomValidity(`Invalid symbol (${mymatch[0]}) at line ${line}, colunm ${iN+1}`);
            return false;
        }
    }
    sequences.setCustomValidity("");
    tokens_span.innerHTML = parseInt(yes_tokens)+'/2048';
    if(yes_tokens < 16){
        sequences.setCustomValidity(`Minimum number of AA allowed: 16. You have ${yes_tokens}`);
        return false;
    }
    if(yes_tokens > 4000){
        sequences.setCustomValidity(`Maximum number of AA allowed: 4000. You have ${yes_tokens}`);
        return false;
    }
    tokens_var.value = yes_tokens;
    return true;
}
window.onload = () => {
    checkValue();
    document.getElementById("sequences").oninput = checkValue;
};
</script>
<style type="text/css" media="screen">
textarea:invalid {
    border: 2px dashed red;
}
textarea:valid {
    border: 2px solid black;
}
input:invalid {
    border: 2px dashed red;
}
input:valid {
    border: 2px solid black;
}
html * { text-align:center; padding:0px !important; }
table, td { border:0px solid black; border-collapse:collapse; }
p { margin:0; }
</style>
</head>
<body>
<h1>Alphafold 2 IBS server</h1>
<form method="post">
    <p>Fasta protein <a href="https://www.uniprot.org/" target="_blank" title="open uniprot for searching">sequence</a>(s):</p>
    <p><textarea id="sequences" name="sequences" rows="20" cols="100" required="true" wrap="off" autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false" style="text-align:left;"></textarea></p>
    <p style="text-align:right"><span>tokens:</span><span id="tokens_span">0/4000</span></p>
    <p>E-mail to get a job-done notification:</p>
    <p><input type="email" id="email" name="email" style="text-align:right;" required autocorrect="off" autocapitalize="off" spellcheck="false" pattern="^[a-z0-9._\-]+@ibs\.fr" placeholder="your_email@ibs.fr"></p>
    <p>&nbsp;</p>
    <input type="submit" value="Submit" style="font-size: 24px">
    <p>&nbsp;</p>
    <input type="hidden" id="tokens_var" name="tokens_var" value="0" />
</form>
<table style="width: 100%">
        <tr>
                <td width="50%">
<p style="text-align:left;">Format accepted:</p>
<p>&nbsp;</p>
<p style="text-align:left;"><b>&gt;descriptor1 (single line)</b></p>
<p style="text-align:left;">FASTA_PROTEIN_SEQUENCE ... (1 or more lines)</p>
<p style="text-align:left;"><b>&gt;descriptor2 (single line)</b></p>
<p style="text-align:left;">FASTA_PROTEIN_SEQUENCE ... (1 or more lines)</p>
<p style="text-align:left;">etc...</p>
<p style="text-align:left;"><b>FOR MULTIMERS</b>:</p>
<p style="text-align:left;"><b>Add as many *copies*</b> of the sequence as you want in the fold</p>
<p style="text-align:left;">WITH <b>*one*</b> &gt;descriptor BEFORE <b>each copy!</b></p>
                </td>
                <td width="50%">
<p style="text-align:left;">Here is a real example:</p>
<pre style="text-align:left;">&gt;WP_267731513.1 AAA family ATPase
MICHAELRSNAKFCDECGGPVAMSSILAEYKQVTV
LFADVVHSMNIAAAVGAERWREIVTELLNRSSTVV
&gt;1LYZ_1|Chain A|HEN EGG WHITE LYSOZYME|Gallus gallus (9031)
KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPC
SALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL
</pre>
                </td>
        </tr>
</table>
    '''
    pagestr += f'''
<h2>Current Slurm job queue:</h2>
<pre>{squeue_output}</pre>
<hr/>
    '''
    pagestr += 'Your IP is: {}\n'.format(client_ip)
    pagestr += '<hr/>\n'
    pagestr += '<h2>Previous Alphafold 2 jobs:</h2>\n'
    pagestr += '<iframe title="Previous AF2 results" width="100%" height="350" src="AF2IBS/previous_results" onload="this.contentWindow.scrollTo(0,this.contentDocument.body.scrollHeight)"></iframe>\n'
    pagestr += '</body>'
    pagestr += '</html>'
    return pagestr

@app.route('/AF2IBS/job/<job_id>', name='job_results')
def job_results(job_id):
    job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)

    if not os.path.exists(job_output_dir):
        return HTTPResponse('Alphafold 2 job not found.', status=404)

    # Check if the job is still running
    finished_flag = os.path.join(job_output_dir, 'input_'+job_id+'_full_dbs/finished.txt')
    running_flag = os.path.join(BASE_OUTPUT_DIR, 'input_'+job_id+'_full_dbs/running.txt')
    failed_flag = os.path.join(job_output_dir, 'input_'+job_id+'_full_dbs/failed.txt')

    if os.path.exists(finished_flag):
        # Job is complete; display results with a link to the file browser
        browse_url = app.get_url('browse', job_id=job_id, filepath='')
        return f'''
            <h1>Alphafold 2 job completed successfully</h1>
            <p>You can browse the output files here:</p>
            <a href="{browse_url}">Browse output files</a>
        '''
    elif os.path.exists(failed_flag):
        # Job has failed; display error message and any available logs
        error_message = ''
        error_log_path = os.path.join(job_output_dir, 'input_'+job_id+'_full_dbs.log')
        if os.path.exists(error_log_path):
            with open(error_log_path, 'r') as f:
                error_message = f.read()
        return f'''
            <h1>Alphafold 2 job failed</h1>
            <p>An error occurred during the job execution.</p>
            <h2>Error Details:</h2>
            <pre>{error_message}</pre>
        '''
    elif os.path.exists(running_flag):
        # Job is still running or did not start yet
        # Get the output of 'nvidia-smi'
        try:
            nvidia_smi_output = subprocess.check_output(['nvidia-smi'], universal_newlines=True)
        except subprocess.CalledProcessError as e:
            nvidia_smi_output = f"Failed to retrieve nvidia-smi output: {e}"
        # Get the output of 'squeue'
        try:
            squeue_output = subprocess.check_output(['squeue'], universal_newlines=True)
        except subprocess.CalledProcessError as e:
            squeue_output = f"Failed to retrieve squeue output: {e}"

        # Return the page indicating the job is still running or in the queue and display nvidia-smi/squeue output
        return f'''
            <h1>Alphafold 2 job is still running or is in the queue.</h1>
            <p>Please refresh this page later.</p>
            <h2>Current Slurm job queue:</h2>
            <pre>{squeue_output}</pre>
            <h2>GPU Status:</h2>
            <pre>{nvidia_smi_output}</pre>
        '''
    else:
        # Job failed to start
        return f'''
            <h1>Alphafold 2 job is not running nor in the queue. (PROBLEM?)</h1>
            <p>Please refresh this page later.</p>
            <p>If this persists, contact leandro.estrozi@ibs.fr</p>
        '''

@app.route('/AF2IBS/browse/<job_id>/<filepath:path>', name='browse_sub')
@app.route('/AF2IBS/browse/<job_id>', name='browse')
def browse(job_id, filepath=''):
    job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)

    if not os.path.exists(job_output_dir):
        return HTTPResponse('Job not found.', status=404)

    # Ensure that the filepath is within the job_output_dir
    safe_root = os.path.abspath(job_output_dir)
    requested_path = os.path.abspath(os.path.join(job_output_dir, filepath))

    if not requested_path.startswith(safe_root):
        return HTTPResponse('Access denied.', status=403)

    if os.path.isdir(requested_path):
        # List directory contents
        items = os.listdir(requested_path)
        items.sort()
        item_links = []
        for item in items:
            item_path = os.path.join(filepath, item)
            fsize = os.path.getsize(requested_path+'/'+item)
#            item_url = app.get_url('browse', job_id=job_id, filepath=item_path)
#            item_url = urllib.parse.quote(item_url)
            if os.path.isdir(os.path.join(requested_path, item)):
                item_url = app.get_url('browse_sub', job_id=job_id, filepath=item_path)
                item_links.append(f'<li><pre style="display:inline">drwxr-xr-x </pre>[<a href="{item_url}">RESULTS</a>] {item}</li>')
            else:
                download_url = app.get_url('download_file', job_id=job_id, filepath=item_path)
                # Check if the file is a PNG image
                if item.lower().endswith('.png'):
                    view_url = app.get_url('view_image', job_id=job_id, filepath=item_path)
                else:
                    view_url = app.get_url('view_file', job_id=job_id, filepath=item_path)
#                    download_url = urllib.parse.quote(download_url)
                item_links.append(f'<li><pre style="display:inline">-rw-r--r-- </pre>[<a href="{download_url}">Download</a>] {item} [<a href="{view_url}" target="_blank">View</a>] {fsize} bytes</li>')

        # Navigation links
        nav_links = ''
        if filepath:
            parent_path = os.path.dirname(filepath)
            parent_url = app.get_url('browse', job_id=job_id, filepath=parent_path)
#            parent_url = urllib.parse.quote(parent_url)
            nav_links = f'<a href="{parent_url}">[Parent Directory]</a>'

        return f'''
            <h1><a href="http://{os.uname()[1]}.ibs.fr:8082/AF2IBS" target="_top">Alphafold 2 IBS server</a></h1>
            <h1>Browsing Alphafold2 results: {filepath or job_id}</h1>
            {nav_links}
            <ul>
                {''.join(item_links)}
            </ul>
        '''
    elif os.path.isfile(requested_path):
        # If a file is requested, serve it for download
        return download_file(job_id, filepath)
    else:
        return HTTPResponse('File or directory not found.', status=404)

@app.route('/AF2IBS/download/<job_id>/<filepath:path>', name='download_file')
def download_file(job_id, filepath):
    job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)

    # Ensure that the filepath is within the job_output_dir
    safe_root = os.path.abspath(job_output_dir)
    requested_file = os.path.abspath(os.path.join(job_output_dir, filepath))

    if not requested_file.startswith(safe_root):
        return HTTPResponse('Access denied.', status=403)

    if not os.path.exists(requested_file):
        return HTTPResponse('File not found.', status=404)

    return static_file(
        os.path.basename(requested_file),
        root=os.path.dirname(requested_file),
        download=os.path.basename(requested_file)
    )

@app.route('/AF2IBS/view_image/<job_id>/<filepath:path>', name='view_image')
def view_image(job_id, filepath):
    job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)

    # Ensure that the filepath is within the job_output_dir
    safe_root = os.path.abspath(job_output_dir)
    requested_file = os.path.abspath(os.path.join(job_output_dir, filepath))

    if not requested_file.startswith(safe_root):
        return HTTPResponse('Access denied.', status=403)

    if not os.path.isfile(requested_file):
        return HTTPResponse('Image not found.', status=404)

    # Serve the image file without download prompt
    return static_file(
        os.path.basename(requested_file),
        root=os.path.dirname(requested_file),
        mimetype='image/png'
    )

@app.route('/AF2IBS/view_file/<job_id>/<filepath:path>', name='view_file')
def view_file(job_id, filepath):
    job_output_dir = os.path.join(BASE_OUTPUT_DIR, job_id)

    # Ensure that the filepath is within the job_output_dir
    safe_root = os.path.abspath(job_output_dir)
    requested_file = os.path.abspath(os.path.join(job_output_dir, filepath))

    if not requested_file.startswith(safe_root):
        return HTTPResponse('Access denied.', status=403)

    if not os.path.isfile(requested_file):
        return HTTPResponse('Image not found.', status=404)

    # Serve the image file without download prompt
    return static_file(
        os.path.basename(requested_file),
        root=os.path.dirname(requested_file),
        mimetype='text/plain'
    )

@app.route('/AF2IBS/')
def wrong():
    redirect("/AF2IBS")

if __name__ == '__main__':
    run(app, host='0.0.0.0', port=8082, debug=True, reloader=True)
#    run(app, host='0.0.0.0', port=8082, debug=False)
