#line 1 "inc/Module/Install/Bundle.pm - /home/ivan/lib/perl5/site_perl/5.8.0/Module/Install/Bundle.pm"
# $File: //depot/cpan/Module-Install/lib/Module/Install/Bundle.pm $ $Author: autrijus $
# $Revision: #7 $ $Change: 1847 $ $DateTime: 2003/12/31 23:14:54 $ vim: expandtab shiftwidth=4

package Module::Install::Bundle;
use Module::Install::Base; @ISA = qw(Module::Install::Base);

use strict;
use Cwd ();
use File::Find ();
use File::Copy ();
use File::Basename ();

sub auto_bundle {
    my $self = shift;

    # Flatten array of arrays into a single array
    my @core = map @$_, map @$_, grep ref, $self->requires;

    $self->bundle(@core);
}

sub bundle {
    my $self = shift;
    $self->admin->bundle(@_) if $self->is_admin;

    my $cwd = Cwd::cwd();
    my $bundles = $self->read_bundles;
    my $bundle_dir = $self->_top->{bundle};
    $bundle_dir =~ s/\W+/\\W+/g;

    while (my ($name, $version) = splice(@_, 0, 2)) {
        next if eval "use $name $version; 1";
        my $source = $bundles->{$name} or next;
        my $target = File::Basename::basename($source);
        mkdir $target or die $! unless -d $target;

        # XXX - clean those directories upon "make clean"?
        File::Find::find({
            wanted => sub {
                my $out = $_;
                $out =~ s/$bundle_dir/./i;
                mkdir $out if -d;
                File::Copy::copy($_ => $out) unless -d;
            },
            no_chdir => 1,
        }, $source);

        $self->bundles($name, $target);
    }

    chdir $cwd;
}

sub read_bundles {
    my $self = shift;
    my %map;

    local *FH;
    open FH, $self->_top->{bundle} . ".yml" or return {};
    while (<FH>) {
        /^(.*?): (['"])?(.*?)\2$/ or next;
        $map{$1} = $3;
    }
    close FH;

    return \%map;
}

sub bundle_deps {
    my ($self, $pkg, $version) = @_;
    my $deps = $self->admin->scan_dependencies($pkg) or return;

    foreach my $key (sort keys %$deps) {
        $self->bundle($key, ($key eq $pkg) ? $version : 0);
    }
}

1;

__END__

#line 172
