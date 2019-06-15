from django.conf.urls import url

from cryptkeeper import views

app_name = 'cryptkeeper'
urlpatterns = [
    url(r'^results/(?P<results_uid>[0-9A-Za-z\-]+)', views.results, name='results'),
    url(r'^', views.get_sequence, name='index'),
]
