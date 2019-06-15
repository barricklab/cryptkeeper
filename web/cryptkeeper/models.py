from django.db import models


import uuid

# Create your models here.
class CryptKeeperResults(models.Model):
    name = models.CharField('Name', max_length=256, blank=True)
    uid = models.UUIDField(default=uuid.uuid1, editable=False, unique=True)
    creation_date = models.DateTimeField('Date', auto_now_add=True)
    #seq_file = FileField('Sequence File', upload_to=owner_directory_path)

    def __str__(self):
        return self.name
