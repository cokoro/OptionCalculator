from django.conf.urls import include, url
#from django.conf.urls import *
from django.contrib import admin
from calc.views import index,add
from calc.g_basket_views import g_basket_index,G_Basket
from calc.a_basket_views import a_basket_index,a_basket
from calc.g_asian_views import g_asian_index,g_asian
from calc.american_views import american_index,american
from calc.european_views import european_index,european
from calc.volatility_views import volatility_index,volatility
admin.autodiscover()

urlpatterns = [
    # Examples:
    url(r'^$', index, name='index'),
    url(r'^a_asian_index/', index, name='index'),
    url(r'^cal/S(.+)K(.+)T(.+)sigma(.+)r(.+)call_put(.+)path(.+)step(.+)controlspecify(.+)/$', add, name='add'),
    url(r'^g_basket_index/', g_basket_index, name='g_basket_index'),
    url(r'^g_basket/S1(.+)S2(.+)sigma1(.+)sigma2(.+)r(.+)T(.+)K(.+)rau(.+)call_put(.+)flag(.+)/$', G_Basket, name='G_Basket'),
    url(r'^a_basket_index/', a_basket_index, name='a_basket_index'),
    url(r'^a_basket/S1(.+)S2(.+)sigma1(.+)sigma2(.+)r(.+)T(.+)K(.+)rau(.+)call_put(.+)flag(.+)/$', a_basket, name='a_basket'),
    url(r'^g_asian_index/', g_asian_index, name='g_asian_index'),
    url(r'^g_asian/S(.+)sigma(.+)r(.+)T(.+)K(.+)step(.+)call_put(.+)/$', g_asian, name='g_asian'),

    url(r'^american_index/', american_index, name='american_index'),
    url(r'^american/S(.+)sigma(.+)r(.+)T(.+)K(.+)N(.+)call_put(.+)/$', american, name='american'),

    url(r'^european_index/', european_index, name='european_index'),
    url(r'^european/S(.+)sigma(.+)r(.+)q(.+)T(.+)K(.+)call_put(.+)/$', european, name='european'),

    url(r'^volatility_index/', volatility_index, name='volatility_index'),
    url(r'^volatility/S(.+)r(.+)q(.+)T(.+)K(.+)premium(.+)call_put(.+)/$', volatility, name='volatility'),

    url(r'^admin/', include(admin.site.urls)),
]
