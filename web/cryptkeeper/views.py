from django.conf import settings
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django import forms
from django.urls import reverse
from django.core.files.base import ContentFile
from django.views.generic import TemplateView

from models import CryptKeeperResults


import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tempfile import NamedTemporaryFile
import os
import subprocess
import uuid
from shutil import copyfile



# Create your views here.
class SeqForm(forms.Form):
    '''
    This defines the DNA sequence submission form
    '''

    sequence_file = forms.FileField(required=False)
    sequence = forms.CharField(widget=forms.Textarea, required=False, min_length=50, max_length=50000)
    organism = forms.ChoiceField(choices=(
        ('ecoli', 'Escherichia coli'),
    #   ('agro', 'Agrobacterium tumefaciens')
    ))

    def clean(self):
        '''
        This is where we detect what type of sequence has been uploaded or placed in in the textbox.
        All form validation happens here.
        '''
        cleaned_data = super(SeqForm, self).clean()
        sequence_file = cleaned_data.get("sequence_file")
        sequence = cleaned_data.get("sequence")

        if sequence_file and sequence:
            raise forms.ValidationError('Please either upload a sequence or copy and paste it into the textbox. '
                                        'Do not submit both.')
        elif not (sequence_file or sequence):
            raise forms.ValidationError('Please submit a sequence for analysis.')
        elif sequence_file:
            if sequence_file.size < 50000:
                with NamedTemporaryFile(delete=False, mode='w') as data_file:
                    for chunk in sequence_file.chunks():
                        data_file.write(chunk)
            else:
                raise forms.ValidationError('File size is too large.')
        elif sequence:
            with NamedTemporaryFile(delete=False, mode='w') as data_file:
                data_file.write(sequence.strip())

        del cleaned_data['sequence']
        del cleaned_data['sequence_file']

        with open(data_file.name, 'rU') as temp_file:
            cleaned_data['temp_file'] = temp_file.name
            lines = temp_file.readlines()
            data = lines[0].strip()

            if re.match(r'^LOCUS', data):
                # Format is genbank...
                try:
                    genome = SeqIO.read(temp_file.name, 'genbank')
                    cleaned_data['features'] = get_genbank_features(genome)
                    cleaned_data['raw_sequence'] = str(genome.seq).upper()
                    cleaned_data['title'] = genome.id
                    output_handle = NamedTemporaryFile(delete=False)
                    SeqIO.write(SeqRecord(Seq(cleaned_data['raw_sequence']),
                                          id=cleaned_data['title']), output_handle, "fasta")
                    output_handle.close()
                    cleaned_data['fasta_file'] = output_handle.name
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: Genbank filetype detected, but file is malformed.')
            elif re.match(r'^\>', data):
                # Format is FASTA...
                try:
                    genome = SeqIO.read(temp_file.name, 'fasta')
                    cleaned_data['raw_sequence'] = str(genome.seq).upper()
                    cleaned_data['title'] = genome.id
                    cleaned_data['fasta_file'] = temp_file.name
                    cleaned_data['features'] = ''
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: FASTA filetype detected, but file is malformed.')

            elif re.match(r'^[ATCGatcg]+', data):
                # Format is just raw sequence
                try:
                    joined_lines = ''.join(lines)
                    cleaned_data['raw_sequence'] = ''.join(joined_lines.split()).upper()
                    cleaned_data['title'] = 'Untitled'
                    cleaned_data['features'] = ''
                    output_handle = NamedTemporaryFile(delete=False)
                    SeqIO.write(SeqRecord(Seq(cleaned_data['raw_sequence']),
                                          id=cleaned_data['title']), output_handle, "fasta")
                    output_handle.close()
                    cleaned_data['fasta_file'] = output_handle.name
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: FASTA file could not be written from textbox input.')

            else:
                if os.path.isfile(cleaned_data['temp_file']):
                    os.remove(cleaned_data['temp_file'])
                raise forms.ValidationError('The submitted sequence is not valid and cannot be processed.')

        return cleaned_data


def get_sequence(request):
    '''
    This view processes the SeqForm defined above.
    '''
    if request.method == 'POST':
        form = SeqForm(request.POST, request.FILES)
        if form.is_valid():

            sequence_file_path = form.cleaned_data.get('fasta_file')
            #ADD: grab name from form

            res = CryptKeeperResults()
            res.name = os.path.basename(sequence_file_path)
            res.save()

            output_path = os.path.join(settings.BASE_DIR, "output")

            if not os.path.exists(output_path):
              os.mkdir(output_path)

            copyfile(sequence_file_path, os.path.join(output_path, str(res.uid) + ".seq"))


            #Experiment with this later
            #Assign file to model
            #f = open('/path/to/hello.world', 'r')
            #res.seq_file = File(f)
            #f.close()

            ck_command = 'cryptkeeper -w -n ' + res.name + ' -i ' + os.path.join(output_path, str(res.uid) + ".seq") + ' -o ' + os.path.join(output_path, str(res.uid))
            print(ck_command)
            subprocess.Popen([ck_command], shell=True)

            os.remove(sequence_file_path)

            #return render(request, 'cryptkeeper/results.html', {'results_name': res.name, 'results_uid': str(res.uid), 'results_div' : None})

            return HttpResponseRedirect(reverse('cryptkeeper:results', kwargs={'results_uid':res.uid} ))

    else:
        form = SeqForm()

    return render(request, 'cryptkeeper/form.html', {'form': form, 'version': '0.1'})



# Show output if it exists
# Otherwise show waiting screen
def results(request, results_uid):

    print("here")
    res = CryptKeeperResults.objects.get(uid=results_uid)

    results_div_path = os.path.join(settings.BASE_DIR, "output", str(res.uid) + ".plot.div")

    results_div_content = None

    if os.path.isfile(results_div_path):
      f = open(results_div_path, 'r')
      results_div_content = f.read()
      f.close()

    return render(request, 'cryptkeeper/results.html', {'results_name': res.name, 'results_uid': res.uid, 'results_div': results_div_content})
